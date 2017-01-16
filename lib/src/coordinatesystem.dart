
part of referencing;

/// The direction of positive increments in the coordinate value for a
/// coordinate system axis. This direction is exact in some cases, and is
/// approximate in other cases.
///
/// Some coordinate systems use non-standard orientations. For example, the
/// first axis in South African grids usually points West, instead of East. This
/// information is obviously relevant for algorithms converting South African
/// grid coordinates into Lat/Long.
enum AxisDirection {
  /// Axis positive direction is towards lower pixel column.
  column_negative,
  /// Axis positive direction is towards higher pixel column.
  column_positive,
  /// Axis positive direction is towards bottom of approximately vertical
  /// display surface.
  display_down,
  /// Axis positive direction is left in display.
  display_left,
  /// Axis positive direction is right in display.
  display_right,
  /// Axis positive direction is towards top of approximately vertical display
  /// surface.
  display_up,
  /// Axis positive direction is down relative to gravity.
  down,
  /// Axis positive direction is π/2 radians clockwise from north.
  east,
  /// Axis positive direction is approximately east-north-east.
  east_north_east,
  /// Axis positive direction is approximately east-south-east.
  east_south_east,
  /// Axis positive direction is towards the future.
  future,
  /// Axis positive direction is in the equatorial plane from the centre of the
  /// modelled earth towards the intersection of the equator with the prime
  /// meridian.
  geocentric_x,
  /// Axis positive direction is in the equatorial plane from the centre of
  /// the modelled earth towards the intersection of the equator and the
  /// meridian π/2 radians eastwards from the prime meridian.
  geocentric_y,
  /// Axis positive direction is from the centre of the modelled earth parallel
  /// to its rotation axis and towards its north pole.
  geocentric_z,
  /// Axis positive direction is north.
  north,
  /// Axis positive direction is approximately north-east.
  north_east,
  /// Axis positive direction is approximately north-north-east.
  north_north_east,
  /// Axis positive direction is approximately north-north-west.
  north_north_west,
  /// Axis positive direction is approximately north-west.
  north_west,
  /// Unknown or unspecified axis orientation.
  other,
  /// Axis positive direction is towards the past.
  past,
  /// Axis positive direction is towards lower pixel row.
  row_negative,
  /// Axis positive direction is towards higher pixel row.
  row_positive,
  /// Axis positive direction is π radians clockwise from north.
  south,
  /// Axis positive direction is approximately south-east.
  south_east,
  /// Axis positive direction is approximately south-south-east.
  south_south_east,
  /// Axis positive direction is approximately south-south-west.
  south_south_west,
  /// Axis positive direction is approximately south-west.
  south_west,
  /// Axis positive direction is up relative to gravity.
  up,
  /// Axis positive direction is 3π/2 radians clockwise from north.
  west,
  /// Axis positive direction is approximately west-north-west.
  west_north_west,
  /// Axis positive direction is approximately west-south-west.
  west_south_west,
}

/// Meaning of the axis value range specified through minimum value and maximum
/// value.
enum RangeMeaning {

  /// Any value between and including minimum value and maximum value is valid.
  exact,

  /// The axis is continuous with values wrapping around at the minimum value
  /// and maximum value.
  wraparound
}

/// Definition of a coordinate system axis.
///
/// This is used to label axes, and indicate the orientation.
class CoordinateSystemAxis {

  /// The abbreviation used for this coordinate system axes.
  final String abbreviation;

  /// Direction of this coordinate system axis.
  final AxisDirection axisDirection;

  /// The unit of measure used for this coordinate system axis.
  final Unit unit;

  /// The minimum value normally allowed for this axis, in the unit of measure
  /// for the axis.
  final double maximumValue;

  /// The maximum value normally allowed for this axis, in the unit of measure
  /// for the axis.
  final double minimumValue;

  /// The meaning of axis value range specified by the minimum and maximum
  /// values.
  final RangeMeaning rangeMeaning;

  /// Constructs an axis.
  ///
  /// When minimum and maximum value are not given, an unbounded axis is
  /// created.
  const CoordinateSystemAxis(this.abbreviation, this.axisDirection, this.unit,
      [this.minimumValue = double.NEGATIVE_INFINITY,
      this.maximumValue = double.INFINITY,
      this.rangeMeaning = RangeMeaning.exact]);


  /// Default axis info for geodetic longitudes in a [GeographicCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are decimal
  /// degrees.
  ///
  /// The ISO 19111 name is `geodetic longitude` and the abbreviation is
  /// `&lambda;` (lambda).
  ///
  /// This axis is usually part of a [GEODETIC_LONGITUDE], [GEODETIC_LATITUDE],
  /// [ELLIPSOIDAL_HEIGHT] set.
  static const CoordinateSystemAxis GEODETIC_LONGITUDE =
  const CoordinateSystemAxis("\u03BB", AxisDirection.east, NonSI.DEGREE_ANGLE);

  /// Default axis info for geodetic latitudes in a [GeographicCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are decimal
  /// degrees.
  ///
  /// The ISO 19111 name is `geodetic latitude` and the abbreviation is `&phi;`
  /// (phi).
  ///
  /// This axis is usually part of a [GEODETIC_LONGITUDE], [GEODETIC_LATITUDE],
  /// [ELLIPSOIDAL_HEIGHT] set.
  static const CoordinateSystemAxis GEODETIC_LATITUDE =
  const CoordinateSystemAxis("\u03C6", AxisDirection.north, NonSI.DEGREE_ANGLE);

  /// Default axis info for longitudes.
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are decimal
  /// degrees.
  ///
  /// The abbreviation is `&lambda;` (lambda).
  ///
  /// This axis is usually part of a [LONGITUDE], [LATITUDE], [ALTITUDE] set.
  static const CoordinateSystemAxis LONGITUDE =
  const CoordinateSystemAxis("\u03BB", AxisDirection.east, NonSI.DEGREE_ANGLE);

  /// Default axis info for latitudes.
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are decimal
  /// degrees.
  ///
  /// The abbreviation is `&phi;` (phi).
  ///
  /// This axis is usually part of a [LONGITUDE], [LATITUDE], [ALTITUDE] set.
  static const CoordinateSystemAxis LATITUDE =
  const CoordinateSystemAxis("\u03C6", AxisDirection.north, NonSI.DEGREE_ANGLE);


  /// The default axis for height values above the ellipsoid in a
  /// [GeographicCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.up] and units are metres.
  ///
  /// The ISO 19111 name is `ellipsoidal height` and the abbreviation is lower
  /// case `h`.
  ///
  /// This axis is usually part of a [GEODETIC_LONGITUDE], [GEODETIC_LATITUDE],
  /// [ELLIPSOIDAL_HEIGHT] set.
  static const CoordinateSystemAxis ELLIPSOIDAL_HEIGHT =
  const CoordinateSystemAxis("h", AxisDirection.up, SI.METER);

  /// The default axis for height values measured from gravity.
  ///
  /// Increasing ordinates values go [AxisDirection.up] and units are metres.
  ///
  /// The ISO 19111 name is `gravity-related height` and the abbreviation is
  /// upper case `H`.
  static const CoordinateSystemAxis GRAVITY_RELATED_HEIGHT =
  const CoordinateSystemAxis("H", AxisDirection.up, SI.METER);

  /// The default axis for altitude values.
  ///
  /// Increasing ordinates values go [AxisDirection.up] and units are metres.
  ///
  /// The abbreviation is lower case `h`.
  ///
  /// This axis is usually part of a [LONGITUDE], [LATITUDE], [ALTITUDE] set.
  static const CoordinateSystemAxis ALTITUDE = const CoordinateSystemAxis(
      "h", AxisDirection.up, SI.METER);

  /// The default axis for depth.
  ///
  /// Increasing ordinates values go [AxisDirection.down] and units are metres.
  ///
  /// The ISO 19111 name is `depth`.
  static const CoordinateSystemAxis DEPTH = const CoordinateSystemAxis(
      "d", AxisDirection.down, SI.METER);

  /// Default axis info for radius in a [GeocentricCRS] using [SphericalCS].
  ///
  /// Increasing ordinates values go [AxisDirection.up] and units are metres.
  ///
  /// The ISO 19111 name is `geocentric radius` and the abbreviation is lower
  /// case `r`.
  ///
  /// This axis is usually part of a [SPHERICAL_LONGITUDE],
  /// [SPHERICAL_LATITUDE], [GEOCENTRIC_RADIUS] set.
  static const CoordinateSystemAxis GEOCENTRIC_RADIUS =
  const CoordinateSystemAxis("r", AxisDirection.up, SI.METER);

  /// Default axis info for longitudes in a [GeocentricCRS] using [SphericalCS].
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are decimal
  /// degrees.
  ///
  /// The ISO 19111 name is `spherical longitude` and the abbreviation is
  /// `&Omega;` (omega).
  ///
  /// This axis is usually part of a [SPHERICAL_LONGITUDE],
  /// [SPHERICAL_LATITUDE], [GEOCENTRIC_RADIUS] set.
  static const CoordinateSystemAxis SPHERICAL_LONGITUDE =
  const CoordinateSystemAxis("\u03A9", AxisDirection.east, NonSI.DEGREE_ANGLE);

  /// Default axis info for latitudes in a [GeocentricCRS] using [SphericalCS].
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are decimal
  /// degrees.
  ///
  /// The ISO 19111 name is `spherical latitude` and the abbreviation is
  /// `&Theta;` (theta).
  ///
  /// This axis is usually part of a [SPHERICAL_LONGITUDE],
  /// [SPHERICAL_LATITUDE], [GEOCENTRIC_RADIUS] set.
  static const CoordinateSystemAxis SPHERICAL_LATITUDE =
  const CoordinateSystemAxis("\u03B8", AxisDirection.north, NonSI.DEGREE_ANGLE);

  /// Default axis info for `x` values in a [CartesianCS].
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are metres.
  ///
  /// The abbreviation is lower case `x`.
  ///
  /// This axis is usually part of a [X], [Y], [Z] set.
  static const CoordinateSystemAxis X = const CoordinateSystemAxis("x",
      AxisDirection.east, SI.METER);

  /// Default axis info for `y` values in a [CartesianCS].
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are metres.
  ///
  /// The abbreviation is lower case `y`.
  ///
  /// This axis is usually part of a [X], [Y], [Z] set.
  static const CoordinateSystemAxis Y = const CoordinateSystemAxis("y",
      AxisDirection.north, SI.METER);

  /// Default axis info for `z` values in a [CartesianCS].
  ///
  /// Increasing ordinates values go [AxisDirection.up] and units are metres.
  ///
  /// The abbreviation is lower case `z`.
  ///
  /// This axis is usually part of a [X], [Y], [Z] set.
  static const CoordinateSystemAxis Z = const CoordinateSystemAxis("z",
      AxisDirection.up, SI.METER);

  /// Default axis info for `x` values in a [GeocentricCRS] using [CartesianCS].
  ///
  /// Increasing ordinates values go toward prime meridian and units are metres.
  ///
  /// The ISO 19111 name is `geocentric X` and the abbreviation is upper case
  /// `X`.
  ///
  /// This axis is usually part of a [GEOCENTRIC_X], [GEOCENTRIC_Y],
  /// [GEOCENTRIC_Z] set.
  static const CoordinateSystemAxis GEOCENTRIC_X = const CoordinateSystemAxis(
      "X", AxisDirection.other, SI.METER);

  /// Default axis info for `y` values in a [GeocentricCRS] using [CartesianCS].
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are metres.
  ///
  /// The ISO 19111 name is `geocentric Y` and the abbreviation is upper case
  /// `Y`.
  ///
  /// This axis is usually part of a [GEOCENTRIC_X], [GEOCENTRIC_Y],
  /// [GEOCENTRIC_Z] set.
  static const CoordinateSystemAxis GEOCENTRIC_Y = const CoordinateSystemAxis(
      "Y", AxisDirection.east, SI.METER);

  /// Default axis info for `z` values in a [GeocentricCRS] using [CartesianCS].
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are metres.
  ///
  /// The ISO 19111 name is `geocentric Z` and the abbreviation is upper case
  /// `Z`.
  ///
  /// This axis is usually part of a [GEOCENTRIC_X], [GEOCENTRIC_Y],
  /// [GEOCENTRIC_Z] set.
  static const CoordinateSystemAxis GEOCENTRIC_Z = const CoordinateSystemAxis(
      "Z", AxisDirection.north, SI.METER);

  /// Default axis info for Easting values in a [ProjectedCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.east] and units are metres.
  ///
  /// The ISO 19111 name is `easting` and the abbreviation is upper case `E`.
  ///
  /// This axis is usually part of a [EASTING], [NORTHING] set.
  static const CoordinateSystemAxis EASTING = const CoordinateSystemAxis("E",
      AxisDirection.east, SI.METER);

  /// Default axis info for Westing values in a [ProjectedCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.west] and units are metres.
  ///
  /// The ISO 19111 name is `westing` and the abbreviation is upper case `W`.
  static const CoordinateSystemAxis WESTING = const CoordinateSystemAxis("W",
      AxisDirection.west, SI.METER);

  /// Default axis info for Northing values in a [ProjectedCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.north] and units are metres.
  ///
  /// The ISO 19111 name is `northing` and the abbreviation is upper case `N`.
  ///
  /// This axis is usually part of a [EASTING], [NORTHING] set.
  static const CoordinateSystemAxis NORTHING = const CoordinateSystemAxis("N",
      AxisDirection.north, SI.METER);

  /// Default axis info for Southing values in a [ProjectedCRS].
  ///
  /// Increasing ordinates values go [AxisDirection.south] and units are metres.
  ///
  /// The ISO 19111 name is `southing` and the abbreviation is upper case `S`.
  static const CoordinateSystemAxis SOUTHING = const CoordinateSystemAxis("S",
      AxisDirection.south, SI.METER);

  /// A default axis for time values in a [TimeCS].
  ///
  /// Increasing time go toward [AxisDirection.future] and units are days.
  ///
  /// The abbreviation is lower case `t`.
  static const CoordinateSystemAxis TIME = const CoordinateSystemAxis("t",
      AxisDirection.future, NonSI.DAY);

  /// A default axis for column indices in a [GridCoverage].
  ///
  /// Increasing values go toward [AxisDirection.column_positive].
  ///
  /// The abbreviation is lower case `i`.
  static const CoordinateSystemAxis COLUMN = const CoordinateSystemAxis("i",
      AxisDirection.column_positive, Unit.one);

  /// A default axis for row indices in a [GridCoverage].
  ///
  /// Increasing values go toward [AxisDirection.row_positive].
  ///
  /// The abbreviation is lower case `j`.
  static const CoordinateSystemAxis ROW = const CoordinateSystemAxis("j",
      AxisDirection.row_positive, Unit.one);

  /// A default axis for `x` values in a display device.
  ///
  /// Increasing values go toward [AxisDirection.display_right].
  ///
  /// The abbreviation is lower case `x`.
  static const CoordinateSystemAxis DISPLAY_X = const CoordinateSystemAxis("x",
      AxisDirection.display_right, Unit.one);

  /// A default axis for `y` values in a display device.
  ///
  /// Increasing values go toward [AxisDirection.display_down].
  ///
  /// The abbreviation is lower case `y`.
  static const CoordinateSystemAxis DISPLAY_Y = const CoordinateSystemAxis("y",
      AxisDirection.display_down, Unit.one);


}

/// The set of coordinate system axes that spans a given coordinate space.
///
/// A coordinate system (CS) is derived from a set of (mathematical) rules for
/// specifying how coordinates in a given space are to be assigned to points.
/// The coordinate values in a coordinate tuple shall be recorded in the order
/// in which the coordinate system axes associations are recorded, whenever
/// those coordinates use a coordinate reference system that uses this
/// coordinate system.
abstract class CoordinateSystem {

  /// The list of axes for this coordinate system.
  final List<CoordinateSystemAxis> axes;

  /// The dimension of the coordinate system.
  int get dimension => axes.length;

  const CoordinateSystem(this.axes);
}


/// A two- or three-dimensional coordinate system in which position is specified
/// by geodetic latitude, geodetic longitude, and (in the three-dimensional
/// case) ellipsoidal height.
///
/// An [EllipsoidalCS] shall have two or three [axes].
class EllipsoidalCS extends CoordinateSystem {

  const EllipsoidalCS(List<CoordinateSystemAxis> axes) : super(axes);

  /// A two-dimensional ellipsoidal CS with
  /// [CoordinateSystemAxis.GEODETIC_LONGITUDE] and
  /// [CoordinateSystemAxis.GEODETIC_LONGITUDE] axis.
  static const EllipsoidalCS GEODETIC_2D = const EllipsoidalCS(const [
    CoordinateSystemAxis.GEODETIC_LONGITUDE,
    CoordinateSystemAxis.GEODETIC_LATITUDE]);

  /// A three-dimensional ellipsoidal CS with
  /// [CoordinateSystemAxis.GEODETIC_LONGITUDE],
  /// [CoordinateSystemAxis.GEODETIC_LONGITUDE] and
  /// [CoordinateSystemAxis.ELLIPSOIDAL_HEIGHT] axis.
  static const EllipsoidalCS GEODETIC_3D = const EllipsoidalCS(const [
    CoordinateSystemAxis.GEODETIC_LONGITUDE,
    CoordinateSystemAxis.GEODETIC_LATITUDE,
    CoordinateSystemAxis.ELLIPSOIDAL_HEIGHT]);


}

/// A two- or three-dimensional coordinate system with straight axes that are
/// not necessarily orthogonal.
///
/// An [AffineCS] shall have two or three axes.
class AffineCS extends CoordinateSystem {
  const AffineCS(List<CoordinateSystemAxis> axes) : super(axes);

}



/// A 1-, 2-, or 3-dimensional coordinate system with straight, orthogonal axes.
///
/// In the multi-dimensional case, all axes shall have the same length unit of
/// measure.
///
/// A [CartesianCS] shall have one, two, or three axes.
class CartesianCS extends AffineCS {
  const CartesianCS(List<CoordinateSystemAxis> axes) : super(axes);


  /// A two-dimensional cartesian CS with Easting and Northing axis in metres.
  static const CartesianCS PROJECTED = const CartesianCS(const [
    CoordinateSystemAxis.EASTING,
    CoordinateSystemAxis.NORTHING]);

  /// A three-dimensional cartesian CS with geocentric x, y and z axis in
  /// metres.
  static const CartesianCS GEOCENTRIC = const CartesianCS(const [
    CoordinateSystemAxis.GEOCENTRIC_X,
    CoordinateSystemAxis.GEOCENTRIC_Y,
    CoordinateSystemAxis.GEOCENTRIC_Z]);

  /// A two-dimensional cartesian CS with x and y axis in metres.
  static const CartesianCS GENERIC_2D = const CartesianCS(const [
    CoordinateSystemAxis.X,
    CoordinateSystemAxis.Y]);

  /// A three-dimensional cartesian CS with x, y and z axis in metres.
  static const CartesianCS GENERIC_3D = const CartesianCS(const [
    CoordinateSystemAxis.X,
    CoordinateSystemAxis.Y,
    CoordinateSystemAxis.Z]);

  /// A two-dimensional cartesian CS with column and row axis.
  static const CartesianCS GRID = const CartesianCS(const [
      CoordinateSystemAxis.COLUMN,
      CoordinateSystemAxis.ROW]);

  /// A two-dimensional cartesian CS with x and y axis.
  static const CartesianCS DISPLAY = const CartesianCS(const [
      CoordinateSystemAxis.DISPLAY_X,
      CoordinateSystemAxis.DISPLAY_Y]);
}

/// A three-dimensional coordinate system consisting of a [PolarCS] extended by
/// a straight coordinate axis perpendicular to the plane spanned by the polar
/// coordinate system.
///
/// A [CylindricalCS] shall have three axes.
class CylindricalCS extends CoordinateSystem {
  const CylindricalCS(List<CoordinateSystemAxis> axes) : super(axes);

}

/// A two-dimensional coordinate system in which position is specified by the
/// distance from the origin and the angle between the line from the origin to a
/// point and a reference direction.
///
/// A [PolarCS] shall have two axes.
class PolarCS extends CoordinateSystem {
  const PolarCS(List<CoordinateSystemAxis> axes) : super(axes);

}

/// A one-dimensional coordinate system that consists of the points that lie on
/// the single axis described.
///
/// The associated ordinate is the distance from the specified origin to the
/// point along the axis. Example: usage of the line feature representing a road
/// to describe points on or along that road. A {@code LinearCS} shall have one
/// [axes].
class LinearCS extends CoordinateSystem {
  const LinearCS(List<CoordinateSystemAxis> axes) : super(axes);

}

/// A three-dimensional coordinate system with one distance measured from the
/// origin and two angular coordinates.
///
/// Not to be confused with an [EllipsoidalCS] based on an ellipsoid
/// "degenerated" into a sphere.
///
/// A [SphericalCS] shall have three axes.
class SphericalCS extends CoordinateSystem {
  const SphericalCS(List<CoordinateSystemAxis> axes) : super(axes);



  /// A three-dimensional spherical CS with
  /// [CoordinateSystemAxis.SPHERICAL_LONGITUDE],
  /// [CoordinateSystemAxis.SPHERICAL_LATITUDE] and
  /// [CoordinateSystemAxis.GEOCENTRIC_RADIUS] axis.
  static const SphericalCS GEOCENTRIC = const SphericalCS(const [
      CoordinateSystemAxis.SPHERICAL_LONGITUDE,
      CoordinateSystemAxis.SPHERICAL_LATITUDE,
      CoordinateSystemAxis.GEOCENTRIC_RADIUS]);


}

/// A one-dimensional coordinate system containing a single time axis, used to
/// describe the temporal position of a point in the specified time units from a
/// specified time origin.
///
/// A [TimeCS] shall have one axis.
class TimeCS extends CoordinateSystem {

  const TimeCS(List<CoordinateSystemAxis> axes) : super(axes);

  /// A one-dimensional temporal CS with time axis in day units.
  static const TimeCS DAYS = const TimeCS(const [CoordinateSystemAxis.TIME]);

  /// A one-dimensional temporal CS with time axis in second units.
  static const TimeCS SECONDS = const TimeCS(const [
    const CoordinateSystemAxis("t", AxisDirection.future, SI.SECOND)
  ]);

}


/// A two- or three-dimensional coordinate system that consists of any
/// combination of coordinate axes not covered by any other Coordinate System
/// type.
///
/// An example is a multilinear coordinate system which contains one coordinate
/// axis that may have any 1-D shape which has no intersections with itself.
/// This non-straight axis is supplemented by one or two straight axes to
/// complete a 2 or 3 dimensional coordinate system. The non-straight axis is
/// typically incrementally straight or curved.
///
/// A [UserDefinedCS] shall have two or three axes.
class UserDefinedCS extends CoordinateSystem {

  const UserDefinedCS(List<CoordinateSystemAxis> axes) : super(axes);

}


/// A one-dimensional coordinate system used to record the heights (or depths)
/// of points.
///
/// Such a coordinate system is usually dependent on the Earth's gravity field,
/// perhaps loosely as when atmospheric pressure is the basis for the vertical
/// coordinate system axis. An exact definition is deliberately not provided as
/// the complexities of the subject fall outside the scope of this
/// specification.
///
/// A [VerticalCS] shall have one axis.
class VerticalCS extends CoordinateSystem {
  const VerticalCS(List<CoordinateSystemAxis> axes) : super(axes);


  /// A one-dimensional vertical CS with
  /// [CoordinateSystemAxis.ELLIPSOIDAL_HEIGHT] axis in metres.
  static const VerticalCS ELLIPSOIDAL_HEIGHT = const VerticalCS(const [
      CoordinateSystemAxis.ELLIPSOIDAL_HEIGHT]);

  /// A one-dimensional vertical CS with
  /// [CoordinateSystemAxis.GRAVITY_RELATED_HEIGHT] axis in metres.
  static const VerticalCS GRAVITY_RELATED_HEIGHT = const VerticalCS(const [
      CoordinateSystemAxis.GRAVITY_RELATED_HEIGHT]);


  /// A one-dimensional vertical CS with [CoordinateSystemAxis.DEPTH] axis in
  /// metres.
  static const VerticalCS DEPTH = const VerticalCS(const [
      CoordinateSystemAxis.DEPTH]);
}