#!/usr/bin/env python3

from math import radians, cos, sin, asin, sqrt, atan2, pi, atan, tan

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon/2) ** 2
    c = 2 * asin(sqrt(a))
    c = c * 180 / pi
    return c

def EarthRadiusInMeters(latitudeRadians):
    a = 6378137.0
    b = 6356752.3
    cos1 = cos(latitudeRadians)
    sin1 = sin(latitudeRadians)
    t1 = a * a * cos1
    t2 = b * b * sin1
    t3 = a * cos1
    t4 = b * sin1
    return sqrt((t1 * t1 + t2 * t2) / (t3 * t3 + t4 * t4))

def GeocentricLatitude(lat):
    e2 = 0.00669437999014
    clat = atan((1.0 - e2) * tan(lat))
    return clat

def LocationToPoint(c, oblate):
    lat = c['lat'] * pi / 180
    lon = c['lon'] * pi / 180
    if (oblate):
        radius = EarthRadiusInMeters(lat)
        clat = GeocentricLatitude(lat)
    else:
        radius = 6371009
        clat = lat

    cosLon = cos(lon)
    sinLon = sin(lon)
    cosLat = cos(clat)
    sinLat = sin(clat)
    x = radius * cosLon * cosLat
    y = radius * sinLon * cosLat
    z = radius * sinLat

    cosGlat = cos(lat)
    sinGlat = sin(lat)

    nx = cosGlat * cosLon
    ny = cosGlat * sinLon
    nz = sinGlat

    x += c['elv'] * nx
    y += c['elv'] * ny
    z += c['elv'] * nz

    return {'x':x, 'y':y, 'z':z, 'radius':radius, 'nx':nx, 'ny':ny, 'nz':nz}

def RotateGlobe(b, a, bradius, aradius, oblate):
    # Get modified coords of b by rotating the globe so that a is at lat=0, lon=0
    br = {'lat':b['lat'], 'lon':(b['lon'] - a['lon']), 'elv':b['elv']}
    brp = LocationToPoint(br, oblate)

    alat = a['lat'] * -1 * pi / 180
    if (oblate):
        alat = GeocentricLatitude(alat)

    acos = cos(alat)
    asin = sin(alat)

    bx = (brp['x'] * acos) - (brp['z'] * asin)
    by = brp['y']
    bz = (brp['x'] * asin) + (brp['z'] * acos)

    return {'x':bx, 'y':by, 'z':bz, 'radius':bradius}

def bearing(lon1, lat1, lon2, lat2):
    """
    Calculate the bearing between two points on the earth
    """

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    a = {'lat':lat1, 'lon':lon1, 'elv':0}
    b = {'lat':lat2, 'lon':lon2, 'elv':0}
    
    ap = LocationToPoint(a, True)
    bp = LocationToPoint(b, True)

    br = RotateGlobe(b, a, bp['radius'], ap['radius'], True)
    if (br['z'] * br['z'] + br['y'] * br['y'] > 0.000001):
        theta = atan2(br['z'], br['y']) * 180.0 / pi
        azimuth = 90.0 - theta
        if (azimuth < 0.0):
            azimuth += 360.0
        if (azimuth > 360.0):
            azimuth -= 360.0

    return theta

#    dlon = lon2 - lon1
    
    # bearing formula
#    b = atan2(sin(dlon) * cos(lat2), cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon))

#    b = b * 180 / pi

#    b = 90 - b
    
#    return b
