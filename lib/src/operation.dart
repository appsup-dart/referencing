
part of referencing;

abstract class CoordinateOperation {

  final CoordinateReferenceSystem sourceCRS;
  @virtual final CoordinateReferenceSystem targetCRS;

  MathTransform get mathTransform;

  const CoordinateOperation(this.sourceCRS, this.targetCRS);

  CoordinateOperation get inverse => new _InverseCoordinateOperation(this);
}

class _CoordinateOperationFromMathTransform extends CoordinateOperation {
  final MathTransform mathTransform;

  const _CoordinateOperationFromMathTransform(
      CoordinateReferenceSystem sourceCRS, CoordinateReferenceSystem targetCRS,
      this.mathTransform) : super(sourceCRS, targetCRS);
}

class _InverseCoordinateOperation extends CoordinateOperation {
  final CoordinateOperation inverse;

  _InverseCoordinateOperation(CoordinateOperation inverse) :
      inverse = inverse, super(inverse.targetCRS, inverse.sourceCRS);

  @override
  MathTransform get mathTransform => inverse.mathTransform.inverse;
}

abstract class SingleOperation extends CoordinateOperation {

  final OperationMethod operationMethod;
  final Map<String,dynamic> parameterValues;

  MathTransform get mathTransform;

  const SingleOperation(CoordinateReferenceSystem sourceCRS,
      CoordinateReferenceSystem targetCRS, this.operationMethod, this.parameterValues) :
      super(sourceCRS, targetCRS);

}


class ConcatenatedOperation extends CoordinateOperation {
  final List<SingleOperation> operations;

  ConcatenatedOperation(List<SingleOperation> operations) :
        operations = operations,
        super(operations.first.sourceCRS, operations.last.targetCRS);

  @override
  MathTransform get mathTransform => new MathTransform.concatenated(
      operations.map((o)=>o.mathTransform));
}



abstract class MathTransform {

  const MathTransform();

  factory MathTransform.concatenated(Iterable<MathTransform> transforms) =>
      transforms.length>1 ? new ConcatenatedTransform(transforms.first,
          new MathTransform.concatenated(transforms.skip(1))) : transforms.first;

  Vector3 transform(Vector3 point);
  Vector3 inverseTransform(Vector3 point);

  MathTransform get inverse => new InverseMathTransform(this);
}

class ConcatenatedTransform extends MathTransform {
  final MathTransform transform1;
  final MathTransform transform2;

  const ConcatenatedTransform(this.transform1, this.transform2);


  @override
  Vector3 inverseTransform(Vector3 point) =>
      transform1.inverseTransform(transform2.inverseTransform(point));

  @override
  Vector3 transform(Vector3 point) =>
      transform2.transform(transform1.transform(point));

}

class HelmertTransformation extends MathTransform {
  final double scaleFactor;

  final double dx;
  final double dy;
  final double dz;

  final double rx;
  final double ry;
  final double rz;

  List<double> get params => [dx,dy,dz,rx,ry,rz,scaleFactor];

  Vector3 get translationVector => new Vector3(dx,dy,dz);
  Matrix3 get rotationMatrix => new Matrix3(1.0, rz, -ry,
      -rz, 1.0, rx,
      ry, -rx, 1.0);
  /*new Matrix3.rotationX(rx)
      * new Matrix3.rotationY(ry)
      * new Matrix3.rotationZ(rz);
  */

  Matrix3 get inverseRotationMatrix => //rotationMatrix.clone()..invert();
  new Matrix3(
      1.0, -rz, ry,
      rz, 1.0, -rx,
      -ry, rx, 1.0
  );
/*
      new Matrix3.rotationZ(-rz)
      * new Matrix3.rotationY(-ry)
      * new Matrix3.rotationX(-rx);
*/

  Matrix4 get matrix => new Matrix4(
      1.0*scaleFactor, -math.sin(rz)*scaleFactor, math.sin(ry)*scaleFactor, dx,
      math.sin(rz)*scaleFactor, 1.0*scaleFactor, -math.sin(rx)*scaleFactor, dy,
      -math.sin(ry)*scaleFactor, math.sin(rx)*scaleFactor, 1.0*scaleFactor, dz,
      0.0, 0.0, 0.0, 1.0);

  const HelmertTransformation(this.dx,this.dy,this.dz,
      double rx,double ry,double rz,double ppm) :
        rx = rx*math.PI/180/3600, ry = ry*math.PI/180/3600, rz = rz*math.PI/180/3600,
        scaleFactor = 1 + ppm/1E+6;

  @override
  Vector3 transform(Vector3 p) =>
      translationVector + rotationMatrix*p*scaleFactor;

  Vector3 inverseTransform(Vector3 p) =>
      inverseRotationMatrix*(p-translationVector)/scaleFactor;


}

class InverseMathTransform extends MathTransform {
  final MathTransform inverse;

  const InverseMathTransform(this.inverse);


  @override
  Vector3 inverseTransform(Vector3 point) => inverse.transform(point);

  @override
  Vector3 transform(Vector3 point) => inverse.inverseTransform(point);

}

abstract class Conversion extends SingleOperation {
  const Conversion(CoordinateReferenceSystem sourceCRS,
      CoordinateReferenceSystem targetCRS, OperationMethod operationMethod, Map<String, dynamic> parameterValues) :
        super(sourceCRS, targetCRS, operationMethod, parameterValues);

}

class Projection extends Conversion {
  const Projection(CoordinateReferenceSystem sourceCRS,
      CoordinateReferenceSystem targetCRS, OperationMethod operationMethod, Map<String, dynamic> parameterValues) :
        super(sourceCRS, targetCRS, operationMethod, parameterValues);

  @override
  ProjectedCRS get targetCRS => super.targetCRS;

  static Expando<MathTransform> _transforms = new Expando();
  @override
  MathTransform get mathTransform => _transforms[this] ??= _createMathTransform();


  MathTransform _createMathTransform() {
    switch (operationMethod.name) {
      case "Lambert_Conformal_Conic_2SP":
        return new LambertConformal.sp2(
            targetCRS.datum.ellipsoid,
            parameterValues["standard_parallel_1"],
            parameterValues["standard_parallel_2"],
            parameterValues["latitude_of_origin"],
            parameterValues["central_meridian"],
            parameterValues["false_easting"],
            parameterValues["false_northing"]
        );
    }
    throw new ArgumentError("Not found");
  }
}

class OperationMethod {

  final String name;

  const OperationMethod(this.name);


}

abstract class MapProjection extends MathTransform {

  final double falseEasting;
  final double falseNorthing;

  const MapProjection(this.falseEasting, this.falseNorthing);
}


class LambertConformal extends MapProjection {

  final Ellipsoid ellipsoid;
  final double standardParallel1;
  final double standardParallel2;
  final double latitudeOfOrigin;
  final double centralMeridian;
  final bool isBelgium;

  double n;
  double f;
  double r_0;

  LambertConformal.sp2(this.ellipsoid, double standardParallel1, double standardParallel2,
      double latitudeOfOrigin, double centralMeridian, double falseEasting,
      double falseNorthing, [this.isBelgium = false]) :
        standardParallel1 = standardParallel1*math.PI/180,
        standardParallel2 = standardParallel2*math.PI/180,
        latitudeOfOrigin = latitudeOfOrigin*math.PI/180,
        centralMeridian = centralMeridian*math.PI/180,
        super(falseEasting, falseNorthing) {
    var lat_1=this.standardParallel1;
    var lat_2=this.standardParallel2;
    var lat_0=this.latitudeOfOrigin;

    var t1 = tsfn(lat_1, math.sin(lat_1));
    var m1 = msfn(math.sin(lat_1), math.cos(lat_1));
    var t2 = tsfn(lat_2, math.sin(lat_2));
    var m2 = msfn(math.sin(lat_2), math.cos(lat_2));

    n = math.log(m1/m2) / math.log(t1/t2);

    f = m1 * math.pow(t1, -n) / n;
//    var f = math.cos(lat_1)*math.pow(math.tan(math.PI/4+lat_1/2),-n)/n; // 11565915.812935
    // f *= ellipsoid.semiMajorAxis;

    r_0 = lat_0==math.PI/2 ? 0.0 :
    f * math.pow(tsfn(latitudeOfOrigin, math.sin(latitudeOfOrigin)), n);//f*math.pow(math.tan(math.PI/4+lat_0/2),-n);
}

  LambertConformal.sp1(Ellipsoid ellipsoid, double latitudeOfOrigin, double centralMeridian,
      double falseEasting, double falseNorthing) : this.sp2(ellipsoid,
      latitudeOfOrigin, latitudeOfOrigin, latitudeOfOrigin, centralMeridian,
      falseEasting, falseNorthing);

  LambertConformal.belgium(Ellipsoid ellipsoid, double standardParallel1,
      double standardParallel2, double latitudeOfOrigin, double centralMeridian,
      double falseEasting, double falseNorthing) : this.sp2(ellipsoid,
      standardParallel1, standardParallel2, latitudeOfOrigin, centralMeridian,
      falseEasting, falseNorthing, true);


  static const BELGE_A = 0.00014204313635987700;

  /**
   * Transforms the specified (<var>&lambda;</var>,<var>&phi;</var>) coordinates
   * (units in radians) and stores the result in {@code ptDst} (linear distance
   * on a unit sphere).
   */
  @override
  Vector3 transform(Vector3 point) {
    var lat = point.y;
    var lon = point.x;
    var lon_0=centralMeridian;

    var r = f*math.pow(tsfn(lat, math.sin(lat)),n);

    var teta = n*(lon-lon_0);
    if (isBelgium) teta -= BELGE_A;

    var globalScale = ellipsoid.semiMajorAxis;

    var x = falseEasting+r*math.sin(teta)*globalScale;
    var y = falseNorthing+(r_0-r*math.cos(teta))*globalScale;

    return new Vector3(x,y,0.0);
  }


  /**
   * Computes function (15-9) and (9-13) from Snyder.
   * Equivalent to negative of function (7-7).
   */
  double tsfn(double phi, double sinphi) {
    var excentricity = math.sqrt(ellipsoid.firstEccentricitySquared);
    sinphi *= excentricity;
    /*
         * NOTE: change sign to get the equivalent of Snyder (7-7).
         */
    return math.tan(0.5 * (math.PI/2 - phi)) / math.pow((1 - sinphi) / (1 + sinphi), 0.5*excentricity);
  }

  /**
   * Computes function <code>f(s,c,e²) = c/sqrt(1 - s²&times;e²)</code> needed for the true scale
   * latitude (Snyder 14-15), where <var>s</var> and <var>c</var> are the sine and cosine of
   * the true scale latitude, and <var>e²</var> is the {@linkplain #excentricitySquared
   * eccentricity squared}.
   */
  double msfn(double s, double c) {
    return c / math.sqrt(1.0 - (s*s) * ellipsoid.firstEccentricitySquared);
  }

  bool get isSpherical => ellipsoid.semiMajorAxis==ellipsoid.semiMinorAxis;
  @override
  Vector3 inverseTransform(Vector3 point) {
    var lon_0=centralMeridian;

    var globalScale = ellipsoid.semiMajorAxis;
    double x = (point.x-falseEasting)/globalScale;
    double y = r_0 - (point.y-falseNorthing)/globalScale;
    double rho = hypot(x, y);  // Zero when the latitude is 90 degrees.
    const EPSILON = 1E-6;
    if (rho > EPSILON) {
      if (n < 0) {
        rho = -rho;
        x = -x;
        y = -y;
      }
      double theta = math.atan2(x, y);
      if (isBelgium) {
        theta += BELGE_A;
      }
      x = theta/n;
      if (isSpherical) {
        y = 2.0 * math.atan(math.pow(f/rho, 1.0/n)) - math.PI/2;
      } else {
        y = cphi2(math.pow(rho/f, 1.0/n));
      }
    } else {
      x = 0.0;
      y = n < 0 ? -(math.PI/2) : (math.PI/2);
    }
    return new Vector3(x+lon_0,y,point.z);
  }

  double hypot(double x,double y)
  {
    double t;
    x = x.abs();
    y = y.abs();
    t = math.min(x,y);
    x = math.max(x,y);
    t = t/x;
    return x*math.sqrt(1+t*t);
  }

  /**
   * Iteratively solve equation (7-9) from Snyder.
   */
  double cphi2(final double ts) {
    const MAXIMUM_ITERATIONS = 15;
    const ITERATION_TOLERANCE = 1E-10;
    var excentricity = math.sqrt(ellipsoid.firstEccentricitySquared);
    final double eccnth = 0.5 * excentricity;
    double phi = (math.PI/2) - 2.0 * math.atan(ts);
    for (int i=0; i<MAXIMUM_ITERATIONS; i++) {
      final double con  = excentricity * math.sin(phi);
      final double dphi = (math.PI/2) - 2.0*math.atan(ts * math.pow((1-con)/(1+con), eccnth)) - phi;
      phi += dphi;
      if (dphi.abs() <= ITERATION_TOLERANCE) {
        return phi;
      }
    }
    throw new Exception("Projection: NO_CONVERGENCE");
  }
}


class GeocentricTransform extends MathTransform {

  final Ellipsoid ellipsoid;

  double get a => ellipsoid.semiMajorAxis;
  double get b => ellipsoid.semiMinorAxis;
  double get e2 => ellipsoid.firstEccentricitySquared;

  const GeocentricTransform(this.ellipsoid);


  @override
  Vector3 inverseTransform(Vector3 point) {
    /* local defintions and variables */
    /* end-criterium of loop, accuracy of sin(Latitude) */
    const genau = 1e-12;
    const genau2 = (genau * genau);
    const maxiter = 30;
    const HALF_PI = math.PI/2;

    var CT; /* sin of geocentric latitude */
    var ST; /* cos of geocentric latitude */
    var RX;
    var RK;
    var RN; /* Earth radius at location */
    var CPHI0; /* cos of start or old geodetic latitude in iterations */
    var SPHI0; /* sin of start or old geodetic latitude in iterations */
    var CPHI; /* cos of searched geodetic latitude */
    var SPHI; /* sin of searched geodetic latitude */
    var SDPHI; /* end-criterium: addition-theorem of sin(Latitude(iter)-Latitude(iter-1)) */
    var iter; /* # of continous iteration, max. 30 is always enough (s.a.) */

    var x = point.x;
    var y = point.y;
    var z = point.z; //Z value not always supplied
    var lon;
    var lat;
    var height;

    var p2 = x * x + y * y;
    var p = math.sqrt(p2); // distance between semi-minor axis and location
    var rr = math.sqrt(p2 + z * z); // distance between center and location


    /*      special cases for latitude and longitude */
    if (p / a < genau) {

      /*  special case, if P=0. (X=0., Y=0.) */
      lon = 0.0;

      /*  if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
     *  of ellipsoid (=center of mass), Latitude becomes PI/2 */
      if (rr / a < genau) {
        lat = HALF_PI;
        height = -b;
        return point;
      }
    } else {
      /*  ellipsoidal (geodetic) longitude
     *  interval: -PI < Longitude <= +PI */
      lon = math.atan2(y, x);
    }

    /* --------------------------------------------------------------
   * Following iterative algorithm was developped by
   * "Institut for Erdmessung", University of Hannover, July 1988.
   * Internet: www.ife.uni-hannover.de
   * Iterative computation of CPHI,SPHI and Height.
   * Iteration of CPHI and SPHI to 10**-12 radian resp.
   * 2*10**-7 arcsec.
   * --------------------------------------------------------------
   */
    CT = z / rr;
    ST = p / rr;
    RX = 1.0 / math.sqrt(1.0 - e2 * (2.0 - e2) * ST * ST);
    CPHI0 = ST * (1.0 - e2) * RX;
    SPHI0 = CT * RX;
    iter = 0;

    /* loop to find sin(Latitude) resp. Latitude
   * until |sin(Latitude(iter)-Latitude(iter-1))| < genau */
    do {
      iter++;
      RN = a / math.sqrt(1.0 - e2 * SPHI0 * SPHI0);

      /*  ellipsoidal (geodetic) height */
      height = p * CPHI0 + z * SPHI0 - RN * (1.0 - e2 * SPHI0 * SPHI0);

      RK = e2 * RN / (RN + height);
      RX = 1.0 / math.sqrt(1.0 - RK * (2.0 - RK) * ST * ST);
      CPHI = ST * (1.0 - RK) * RX;
      SPHI = CT * RX;
      SDPHI = SPHI * CPHI0 - CPHI * SPHI0;
      CPHI0 = CPHI;
      SPHI0 = SPHI;
    }
    while (SDPHI * SDPHI > genau2 && iter < maxiter);

    /*      ellipsoidal (geodetic) latitude */
    lat = math.atan(SPHI / (CPHI).abs());
    return new Vector3(lon, lat, height);
  }

  /**
   * Converts geodetic coordinates (longitude, latitude, height) to geocentric
   * coordinates (x, y, z) according to the current ellipsoid parameters.
   */
  @override
  Vector3 transform(Vector3 point) {
    var lon = point.x;
    var lat = point.y;
    var height = point.z;

    const HALF_PI = math.PI/2;
    /*
   ** Don't blow up if Latitude is just a little out of the value
   ** range as it may just be a rounding issue.  Also removed longitude
   ** test, it should be wrapped by Math.cos() and Math.sin().  NFW for PROJ.4, Sep/2001.
   */
    if (lat < -HALF_PI && lat > -1.001 * HALF_PI) {
      lat = -HALF_PI;
    } else if (lat > HALF_PI && lat < 1.001 * HALF_PI) {
      lat = HALF_PI;
    } else if ((lat < -HALF_PI) || (lat > HALF_PI)) {
      throw new ArgumentError("Latitude out of range");
    }

    if (lon > math.PI) {
      lon -= (2 * math.PI);
    }


    var sinLat = math.sin(lat);
    var cosLat = math.cos(lat);
    var sinLat2 = sinLat * sinLat;
    var Rn = a / (math.sqrt(1.0e0 - e2 * sinLat2));
    return new Vector3(
        (Rn + height) * cosLat * math.cos(lon),
        (Rn + height) * cosLat * math.sin(lon),
        ((Rn * (1 - e2)) + height) * sinLat
    );
  }
}

class AngleUnitTransform extends MathTransform {

  final UnitConverter converter;

  const AngleUnitTransform(this.converter);

  @override
  Vector3 inverseTransform(Vector3 point) {
    var converter = this.converter.inverse;
    return new Vector3(
        converter.convert(point.x),
        converter.convert(point.y),
        point.z
    );
  }

  @override
  Vector3 transform(Vector3 point) {
    return new Vector3(
        converter.convert(point.x),
        converter.convert(point.y),
        point.z
    );
  }
}