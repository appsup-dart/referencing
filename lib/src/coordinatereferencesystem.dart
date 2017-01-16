
part of referencing;


abstract class ReferenceSystem {
  const ReferenceSystem();
}

abstract class CoordinateReferenceSystem extends ReferenceSystem {

  @virtual final CoordinateSystem coordinateSystem;

  const CoordinateReferenceSystem(this.coordinateSystem);


  CoordinateOperation conversionTo(CoordinateReferenceSystem target) {
    if (this==target) {
      // TODO
    }
    if (target is ProjectedCRS) {
      return new ConcatenatedOperation(<SingleOperation>[
        conversionTo(target.baseCRS),
        target.conversionFromBase
      ]);
    }
    if (this is ProjectedCRS) {
      return new ConcatenatedOperation(<SingleOperation>[
        (this as ProjectedCRS).conversionFromBase.inverse,
        (this as ProjectedCRS).baseCRS.conversionTo(target)]);
    }
    if (this is GeographicCRS&&target is GeographicCRS) {
      return new _CoordinateOperationFromMathTransform(
          this, target, new MathTransform.concatenated([
        new AngleUnitTransform((this as GeographicCRS).datum.primeMeridian.angularUnit.getConverterTo(SI.RADIAN)),
        new GeocentricTransform((this as GeographicCRS).datum.ellipsoid),
        (this as GeographicCRS).datum.toWgs84Helmert,
        target.datum.toWgs84Helmert.inverse,
        new GeocentricTransform(target.datum.ellipsoid).inverse,
        new AngleUnitTransform(target.datum.primeMeridian.angularUnit.getConverterTo(SI.RADIAN).inverse),
      ]));
    }
    throw new UnimplementedError();
  }
}

class GeodeticCRS extends SingleCRS {

  const GeodeticCRS(GeodeticDatum datum, CoordinateSystem coordinateSystem) :
        super(datum, coordinateSystem);

}

class GeographicCRS extends GeodeticCRS {

  static const GeographicCRS WGS84 = const GeographicCRS(
      GeodeticDatum.WGS84, EllipsoidalCS.GEODETIC_2D);


  static const GeographicCRS WGS84_3D = const GeographicCRS(
      GeodeticDatum.WGS84, EllipsoidalCS.GEODETIC_3D);

  static const GeographicCRS BELGE_1972 = const GeographicCRS(
      GeodeticDatum.BELGE_1972, EllipsoidalCS.GEODETIC_3D);

  @override
  EllipsoidalCS get coordinateSystem => super.coordinateSystem;

  const GeographicCRS(GeodeticDatum datum, EllipsoidalCS coordinateSystem) :
        super(datum, coordinateSystem);
}

class GeocentricCRS extends GeodeticCRS {


  static const GeocentricCRS CARTESIAN = const GeocentricCRS(
      GeodeticDatum.WGS84, CartesianCS.GEOCENTRIC);


  const GeocentricCRS(GeodeticDatum datum, CoordinateSystem coordinateSystem) :
        super(datum, coordinateSystem);
}

/// Abstract coordinate reference system, consisting of a single Coordinate
/// System and a single Datum (as opposed to Compound CRS).
///
/// A coordinate reference system consists of an ordered sequence of coordinate
/// system axes that are related to the earth through a datum. A coordinate
/// reference system is defined by one datum and by one coordinate system. Most
/// coordinate reference system do not move relative to the earth, except for
/// engineering coordinate reference systems defined on moving platforms such
/// as cars, ships, aircraft, and spacecraft.
///
/// Coordinate reference systems are commonly divided into sub-types. The
/// common classification criterion for sub-typing of coordinate reference
/// systems is the way in which they deal with earth curvature. This has a
/// direct effect on the portion of the earth's surface that can be covered by
/// that type of CRS with an acceptable degree of error. The exception to the
/// rule is the subtype "Temporal" which has been added by analogy.
class SingleCRS extends CoordinateReferenceSystem {
  final GeodeticDatum datum;

  const SingleCRS(this.datum, CoordinateSystem coordinateSystem) : super(coordinateSystem);
}


abstract class GeneralDerivedCRS extends SingleCRS {

  final CoordinateReferenceSystem baseCRS;

  Conversion get conversionFromBase;

  const GeneralDerivedCRS(this.baseCRS,
      GeodeticDatum datum, CoordinateSystem coordinateSystem) :
        super(datum, coordinateSystem);

}

class ProjectedCRS extends GeneralDerivedCRS {

  Projection get conversionFromBase => new Projection(baseCRS, this, _operationMethod, _parameterValues);

  final OperationMethod _operationMethod;
  final Map<String,dynamic> _parameterValues;

  const ProjectedCRS(CoordinateReferenceSystem baseCRS,
      this._operationMethod, this._parameterValues, GeodeticDatum datum,
      CoordinateSystem coordinateSystem) :
        super(baseCRS, datum, coordinateSystem);

  /*
PROJCS["Belge 1972 / Belgian Lambert 72",
    GEOGCS["Belge 1972",
        DATUM["Reseau_National_Belge_1972",
            SPHEROID["International 1924",6378388,297,
                AUTHORITY["EPSG","7022"]],
            TOWGS84[-106.869,52.2978,-103.724,0.3366,-0.457,1.8422,-1.2747],
            AUTHORITY["EPSG","6313"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4313"]],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    PROJECTION["Lambert_Conformal_Conic_2SP"],
    PARAMETER["standard_parallel_1",51.16666723333333],
    PARAMETER["standard_parallel_2",49.8333339],
    PARAMETER["latitude_of_origin",90],
    PARAMETER["central_meridian",4.367486666666666],
    PARAMETER["false_easting",150000.013],
    PARAMETER["false_northing",5400088.438],
    AUTHORITY["EPSG","31370"],
    AXIS["X",EAST],
    AXIS["Y",NORTH]]
 */
  static const ProjectedCRS BELGE_1972 = const ProjectedCRS(
      GeographicCRS.BELGE_1972,
      const OperationMethod("Lambert_Conformal_Conic_2SP"),
      const {
        "standard_parallel_1":51.16666723333333,
        "standard_parallel_2":49.8333339,
        "latitude_of_origin":90.0,
        "central_meridian":4.367486666666666,
        "false_easting":150000.013,
        "false_northing":5400088.438,
      },
      GeodeticDatum.BELGE_1972,
      CartesianCS.PROJECTED
  );
}
