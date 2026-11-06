"""Minimal stand-ins for the Orekit/JPype modules Cell 1 imports, so its pure-
Python parts (ROI definitions, WDPA loading, H3 tiling, summaries) can run
here, where orekit-data cannot be downloaded. Nothing here is used for any
orbit or frame computation in the tests."""
import sys, types, math


def install():
    def mod(name, **attrs):
        m = types.ModuleType(name)
        m.__dict__.update(attrs)
        sys.modules[name] = m
        return m

    class FastMath:
        toRadians = staticmethod(math.radians)
        toDegrees = staticmethod(math.degrees)

    class Constants:
        WGS84_EARTH_EQUATORIAL_RADIUS = 6378137.0
        WGS84_EARTH_FLATTENING = 1.0 / 298.257223563

    class IERSConventions:
        IERS_2010 = "IERS_2010"

    class FramesFactory:
        @staticmethod
        def getITRF(*a): return "ITRF(stub)"

    class OneAxisEllipsoid:
        def __init__(self, *a): self.args = a

    class GeodeticPoint:
        def __init__(self, lat, lon, alt): self.lat, self.lon, self.alt = lat, lon, alt

    class TopocentricFrame:
        def __init__(self, *a): pass

    mod("jdk4py", JAVA_HOME="/nonexistent")
    mod("orekit_jpype", initVM=lambda *a, **k: None)
    mod("orekit_jpype.pyhelpers", setup_orekit_curdir=lambda **k: None, setup_orekit_data=lambda **k: None)
    for p in ("org", "org.orekit", "org.hipparchus"):
        mod(p)
    mod("org.orekit.frames", FramesFactory=FramesFactory, TopocentricFrame=TopocentricFrame)
    mod("org.orekit.utils", IERSConventions=IERSConventions, Constants=Constants)
    mod("org.orekit.bodies", OneAxisEllipsoid=OneAxisEllipsoid, GeodeticPoint=GeodeticPoint)
    mod("org.hipparchus.util", FastMath=FastMath)
