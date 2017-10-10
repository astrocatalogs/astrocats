"""
"""
import os
import datetime

import filecmp
import urllib
from bs4 import BeautifulSoup

from astropy import units as un
from astropy.coordinates import SkyCoord

from . import PATH_SDSS_MISSING_HOST_IMAGE, PATH_SDSS_FAILED_HOST_IMAGE
from . producer import Producer_Base

# Force all host images to be reloaded
FORCE_RELOAD_IMAGES_FLAG = False

_IMAGE_SIZE_REQUEST = 512

# Skyservice is used to get SDSS images
# -------------------------------------
# Options:  G: Grid, L: Label, P: PhotoObjs, S: SpecObjs, T: TargetObjs, O: Outline,
#           B: BoundingBox, F: Fields, M: Masks, Q: Plates, I: InvertImage
_SKYSERVICE_OPTIONAL_FLAGS = "GL"

# See: http://skyserver.sdss.org/dr12/en/help/docs/api.aspx#imgcutout
# Example: ("http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg?"
#           "ra=224.5941&dec=-1.09&width=512&height=512&opt=OG")
_URL_SKYSERVICE_IMAGE_REQUEST = (
    # "http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx?"  # NOTE: May2017 not working
    "http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg?"
    "ra={RA_DEG_STR}&dec={DEC_DEG_STR}&scale={SCALE}&width={WIDTH}&height={HEIGHT}&opt={FLAGS}")

_URL_SKYSERVICE_IMAGE_LINK = ("http://skyserver.sdss.org/DR13/en/tools/chart/navi.aspx?opt=G&"
                              "ra={RA_DEG_STR}&dec={DEC_DEG_STR}&scale=0.15")


# Skyview is used to get DSS images
# ---------------------------------
_URL_SKYVIEW_IMAGE_REQUEST = (
    "http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Position={RA_DEC}&coordinates=J2000"
    "&coordinates=&projection=Tan&pixels={SIZE}&size={SCALE}&float=on&scaling=Log&"
    "resolver=SIMBAD-NED&Sampler=_skip_&Deedger=_skip_&rotation=&Smooth=&"
    "lut=colortables%2Fb-w-linear.bin&PlotColor=&grid=_skip_&gridlabels=1&"
    "catalogurl=&CatalogIDs=on&RGB=1&survey=DSS2+IR&survey=DSS2+Red&survey=DSS2+Blue&"
    "IOSmooth=&contour=&contourSmooth=&ebins=null"
)


class HOST_IMAGE_KEYS:
    # Required information for constructing webpages
    IMAGE_PATH = "image_path"
    LINK_URL = "link_url"

    # Auxilliary information
    SOURCE = "source"
    QUERY_URL = "query_url"
    DOWNLOAD_DATE = "download_date"


class Host_Image_Pro(Producer_Base):
    """
    """

    SDSS_IMAGE_SCALE = 0.3
    DSS_IMAGE_SCALE = 0.13889 * SDSS_IMAGE_SCALE

    RUN_FINISH_ON_CHECKPOINT = True

    def __init__(self, catalog, args):
        log = catalog.log
        self.log = log
        log.debug("Host_Image_Pro.__init__()")

        self.host_image_fname = catalog.PATHS.HOST_IMAGES_FILE
        self.dir_images = catalog.PATHS.PATH_HTML
        log.debug("`host_image_fname` = '{}'".format(self.host_image_fname))
        log.debug("`dir_images` = '{}'".format(self.dir_images))

        self.host_images = self._load_default_json(self.host_image_fname)
        log.info("Loaded {} host image locations.".format(len(self.host_images)))
        self.count_exist = 0
        self.count_added = 0
        return

    def update(self, fname, event_name, event_data):
        """Retrieve and/or update the host-images dictionary for this event.

        Dictionary entries have keys and values from `HOST_IMAGE_KEYS`.

        Description
        -----------
        If an entry for this event exists, load it.
        If not, try to download host images from online.
            A new entry is created regardless.
            If images are downloaded, `HOST_IMAGE_KEYS` parameters are added to the new entry.
            If images are *not* downloaded, then only the datetime is recorded.
        If images exist, then their file-location and link-url are returned (used by HTML_Pro).

        If the event does *not* have RA/Dec, or they are erroneous, then no entry is added.
        Otherwise an entry should always exist if the event has been processes.

        Returns
        -------
        None if no images exist
        (event_image_path, link_url) if images are found

        """
        self.log.debug("Host_Image_Pro.update()")

        if ('ra' not in event_data) or ('dec' not in event_data):
            self.log.info("event '{}' missing ra/dec, skipping host image.".format(event_name))
            return

        # Convert / Process Coordinates
        # -----------------------------
        ra_event = event_data['ra'][0]['value']
        dec_event = event_data['dec'][0]['value']
        try:
            coord = SkyCoord(ra=ra_event, dec=dec_event, unit=(un.hourangle, un.deg))
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception:
            warn_str = "Malformed angle for event '{}' - ra: '{}', dec: '{}'".format(
                event_name, ra_event, dec_event)
            self.log.error(warn_str)
            return

        now = datetime.datetime.now().strftime("%Y/%m/%d - %H:%M:%S")

        # Load existing entry for this event or construct a new one
        # ---------------------------------------------------------
        if (event_name in self.host_images) and (not FORCE_RELOAD_IMAGES_FLAG):
            host_image_entry = self.host_images[event_name]
        else:

            # Try to get SDSS Image
            host_image_entry = self._download_sdss_image(event_name, coord)

            # Get DSS Image (on SDSS failure)
            if host_image_entry is None:
                host_image_entry = self._download_dss_image(event_name, ra_event, dec_event)

            if host_image_entry is None:
                host_image_entry = {}
            else:
                self.count_added += 1

            # Store the time of download
            host_image_entry[HOST_IMAGE_KEYS.DOWNLOAD_DATE] = now

            # Store dictionary entry for this event
            self.host_images[event_name] = host_image_entry

        # If we dont have the image, return None
        if (((HOST_IMAGE_KEYS.LINK_URL not in host_image_entry) or
             (HOST_IMAGE_KEYS.IMAGE_PATH not in host_image_entry))):
            return

        # Extract the required parameters from the dictionary entry
        # ---------------------------------------------------------
        link_url = host_image_entry[HOST_IMAGE_KEYS.LINK_URL]
        event_image_path = host_image_entry[HOST_IMAGE_KEYS.IMAGE_PATH]
        # Make sure the stored image path exists
        if not os.path.exists(event_image_path):
            err = "Event '{}': path '{}' does not exist!".format(event_name, event_image_path)
            self.log.error(err)
            # Delete erroneous entry
            del self.host_images[event_name]
            return

        self.count_exist += 1

        return event_image_path, link_url

    def finish(self, lvl=None):
        self.log.info("{} Images exist.  {} Added.".format(self.count_exist, self.count_added))
        self._save_json(self.host_image_fname, self.host_images, expanded=True, lvl=lvl)
        return

    def _get_image_fname(self, event_name, source):
        fname = "{EVENT_NAME}__host_{SOURCE}.jpg".format(EVENT_NAME=event_name, SOURCE=source)
        path = os.path.join(self.dir_images, fname)
        return fname, path

    def _download_sdss_image(self, event_name, coord):
        """Try to download the host for this event from SDSS (via: 'skyservice')

        Returns
        -------
        host_image_entry : dict or None
            On failure, `None` is returned.
            On success, a dict with keys matching `HOST_IMAGE_KEYS` is returned

        """
        ra_str = str(coord.ra.deg)
        dec_str = str(coord.dec.deg)
        SRC = "SDSS-DR13"
        event_image_fname, event_image_path = self._get_image_fname(event_name, SRC)

        # Construct the query URL
        query_url = _URL_SKYSERVICE_IMAGE_REQUEST.format(
            RA_DEG_STR=ra_str, DEC_DEG_STR=dec_str, SCALE=self.SDSS_IMAGE_SCALE,
            WIDTH=_IMAGE_SIZE_REQUEST, HEIGHT=_IMAGE_SIZE_REQUEST,
            FLAGS=_SKYSERVICE_OPTIONAL_FLAGS)
        # Query server
        self.log.debug("query_url: '{}'".format(query_url))

        try:
            response = urllib.request.urlopen(query_url, timeout=60)
            resptxt = response.read()
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as err:
            self.log.info("Event '{}' query '{}' failed: '{}'".format(
                event_name, query_url, str(err)))
            return

        # Save image
        self._save(event_image_path, resptxt, format='wb', lvl=self.log.INFO)

        # Compare image to 'missing' and 'failed' host-images (retrieved if no image exists)
        if filecmp.cmp(event_image_path, PATH_SDSS_MISSING_HOST_IMAGE):
            self.log.info("Missing SDSS image for '{}'".format(event_name))
            return

        if filecmp.cmp(event_image_path, PATH_SDSS_FAILED_HOST_IMAGE):
            self.log.info("Failed SDSS image for '{}'".format(event_name))
            return

        # Construct the URL to link to skyservice (NOTE: different than query url)
        link_url = _URL_SKYSERVICE_IMAGE_LINK.format(RA_DEG_STR=ra_str, DEC_DEG_STR=dec_str)

        host_image_entry = {
            HOST_IMAGE_KEYS.SOURCE: SRC,
            HOST_IMAGE_KEYS.QUERY_URL: query_url,
            HOST_IMAGE_KEYS.LINK_URL: link_url,
            HOST_IMAGE_KEYS.IMAGE_PATH: event_image_path,
        }

        return host_image_entry

    def _download_dss_image(self, event_name, ra_event, dec_event):
        """
        """
        SRC = "DSS-2"
        host_image_entry = None
        event_image_fname, event_image_path = self._get_image_fname(event_name, SRC)

        try:
            ra_dec = str(urllib.parse.quote_plus(ra_event + " " + dec_event))
            temp_url = _URL_SKYVIEW_IMAGE_REQUEST.format(
                RA_DEC=ra_dec, SCALE=self.DSS_IMAGE_SCALE, SIZE=_IMAGE_SIZE_REQUEST)
            self.log.debug("skyview_url: '{}'".format(temp_url))
            response = urllib.request.urlopen(temp_url, timeout=60)
            bandsoup = BeautifulSoup(response, "html5lib")
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as err:
            self.log.warning("Event '{}' query '{}' failed: '{}'".format(
                event_name, temp_url, str(err)))
            return host_image_entry

        images = bandsoup.findAll('img')
        imgname = ''
        for image in images:
            if "Quicklook RGB image" in image.get('alt', ''):
                imgname = image.get('src', '').split('/')[-1]
                break

        if imgname:
            query_url = 'http://skyview.gsfc.nasa.gov/tempspace/fits/' + imgname
            try:
                response = urllib.request.urlopen(query_url)
                self._save(event_image_path, response.read(), format='wb', lvl=self.log.INFO)
                self.count_added += 1
                # local_url = urllib.parse.quote(event_image_fname)
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as err:
                self.log.warning("Event '{}' query '{}' failed: '{}'".format(
                    event_name, query_url, str(err)))
                return host_image_entry

            host_image_entry = {
                HOST_IMAGE_KEYS.SOURCE: SRC,
                HOST_IMAGE_KEYS.QUERY_URL: query_url,
                HOST_IMAGE_KEYS.LINK_URL: query_url,
                HOST_IMAGE_KEYS.IMAGE_PATH: event_image_path,
            }

        return host_image_entry
