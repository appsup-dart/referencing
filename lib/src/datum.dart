
part of referencing;


class Datum {

  const Datum();

}

class GeodeticDatum extends Datum {

  static const GeodeticDatum WGS84 = const GeodeticDatum(
      Ellipsoid.WGS84, PrimeMeridian.GREENWICH,
      const HelmertTransformation(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

  static const GeodeticDatum BELGE_1972 = const GeodeticDatum(
      Ellipsoid.INTERNATIONAL_1924, PrimeMeridian.GREENWICH_RADIANS,
      const HelmertTransformation(-106.868628,52.297783,-103.723893,0.336570,
          -0.456955,1.842183,-1.2747));

  final Ellipsoid ellipsoid;
  final PrimeMeridian primeMeridian;

  final HelmertTransformation toWgs84Helmert;

  const GeodeticDatum(this.ellipsoid, this.primeMeridian, this.toWgs84Helmert);

}


class Ellipsoid {

  final double semiMajorAxis;
  final double inverseFlattening;

  final double flattening;
  final double semiMinorAxis;
  final double firstEccentricitySquared;
  final double secondEccentricitySquared;

  const Ellipsoid(double semiMajorAxis, double inverseFlattening) :
      this._fromSemiMajorAxisAndFlattening(semiMajorAxis, 1/inverseFlattening);

  const Ellipsoid._fromSemiMajorAxisAndFlattening(double semiMajorAxis, double flattening) :
        this.semiMajorAxis = semiMajorAxis,
        this.inverseFlattening = 1/flattening,
        flattening = flattening,
        semiMinorAxis = semiMajorAxis*(1-flattening),
        firstEccentricitySquared = flattening*(2-flattening),
        secondEccentricitySquared = flattening*(2-flattening)/(1-flattening)/(1-flattening);


  static const Ellipsoid WGS84 = const Ellipsoid(6378137.0, 298.257223563);
  static const Ellipsoid INTERNATIONAL_1924 = const Ellipsoid(6378388.0, 297.0);

}

class PrimeMeridian {
  final Unit<Angle> angularUnit;

  final greenwichLongitude;

  const PrimeMeridian(this.greenwichLongitude, [this.angularUnit = NonSI.DEGREE_ANGLE]);

  static const PrimeMeridian GREENWICH = const PrimeMeridian(0.0);
  static const PrimeMeridian GREENWICH_RADIANS = const PrimeMeridian(0.0, SI.RADIAN);
}

