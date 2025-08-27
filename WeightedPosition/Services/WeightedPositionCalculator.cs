//using Mapsui.Projections;
//using Mapsui;
//using NetTopologySuite.Geometries;
//using WeightedPosition.Models;

//namespace WeightedPosition.Services
//{
//    public class WeightedPositionCalculator
//    {
//        private readonly GeometryFactory _factory = new();

//        public MPoint? CalculateWeightedCenter(IEnumerable<Station> stations)
//        {
//            var stationList = stations.Where(s => s.Quality > 0).ToList();
//            if (stationList.Count < 2) return null;

//            // 1. Проецируем станции
//            var projected = stationList.Select(s => new
//            {
//                Station = s,
//                Point = Project(s.LongitudeDeg, s.LatitudeDeg)
//            }).ToList();

//            // 2. Строим лучи
//            var rays = projected.Select(p => new
//            {
//                Ray = BuildRay(p.Point, p.Station.BearingDeg),
//                Weight = p.Station.Quality
//            }).ToList();

//            // 3. Находим пересечения
//            var intersections = new List<(Point point, double weight)>();
//            for (int i = 0; i < rays.Count; i++)
//            {
//                for (int j = i + 1; j < rays.Count; j++)
//                {
//                    var inter = rays[i].Ray.Intersection(rays[j].Ray);
//                    if (inter is Point p)
//                    {
//                        double avgWeight = (rays[i].Weight + rays[j].Weight) / 2.0;
//                        intersections.Add((p, avgWeight));
//                    }
//                }
//            }

//            if (intersections.Count == 0) return null;

//            // 4. Считаем взвешенный центр
//            double sumX = 0, sumY = 0, sumW = 0;
//            foreach (var (pt, w) in intersections)
//            {
//                sumX += pt.X * w;
//                sumY += pt.Y * w;
//                sumW += w;
//            }

//            if (sumW == 0) return null;

//            return new MPoint(sumX / sumW, sumY / sumW);
//        }

//        private MPoint Project(double lon, double lat)
//        {
//            var (x, y) = SphericalMercator.FromLonLat(lon, lat);
//            return new MPoint(x, y);
//        }

//        private LineString BuildRay(MPoint origin, double bearingDeg, double length = 200_000)
//        {
//            double angleRad = (90 - bearingDeg) * Math.PI / 180.0;
//            var dx = length * Math.Cos(angleRad);
//            var dy = length * Math.Sin(angleRad);

//            var start = new Coordinate(origin.X, origin.Y);
//            var end = new Coordinate(origin.X + dx, origin.Y + dy);
//            return new LineString(new[] { start, end });
//        }
//    }
//}
using Mapsui;
using Mapsui.Projections;
using NetTopologySuite.Geometries;
using WeightedPosition.Models;

namespace WeightedPosition.Services
{
    public class WeightedPositionCalculator
    {
        private readonly GeometryFactory _factory = new();

        public MPoint? CalculateWeightedCenter(IEnumerable<Station> stations)
        {
            var stationList = stations.Where(s => s.Quality > 0).ToList();
            if (stationList.Count < 2) return null;

            // 1. Проецируем станции в Меркатор для построения лучей
            var projected = stationList.Select(s => new
            {
                Station = s,
                Point = Project(s.LongitudeDeg, s.LatitudeDeg)
            }).ToList();

            // 2. Строим лучи
            var rays = projected.Select(p => new
            {
                Ray = BuildRay(p.Point, p.Station.BearingDeg),
                Weight = p.Station.Quality
            }).ToList();

            // 3. Находим пересечения
            var intersections = new List<(double lon, double lat, double weight)>();
            for (int i = 0; i < rays.Count; i++)
            {
                for (int j = i + 1; j < rays.Count; j++)
                {
                    var inter = rays[i].Ray.Intersection(rays[j].Ray);
                    if (inter is Point p)
                    {
                        double avgWeight = (rays[i].Weight + rays[j].Weight) / 2.0;

                        // Переводим пересечение из Меркатора обратно в WGS84
                        var (lon, lat) = SphericalMercator.ToLonLat(p.X, p.Y);

                        intersections.Add((lon, lat, avgWeight));
                    }
                }
            }

            if (intersections.Count == 0) return null;

            // 4. Считаем взвешенный центр на сфере через 3D-вектора
            double sumX = 0, sumY = 0, sumZ = 0, sumW = 0;
            foreach (var (lon, lat, w) in intersections)
            {
                double latRad = DegToRad(lat);
                double lonRad = DegToRad(lon);

                double x = Math.Cos(latRad) * Math.Cos(lonRad);
                double y = Math.Cos(latRad) * Math.Sin(lonRad);
                double z = Math.Sin(latRad);

                sumX += x * w;
                sumY += y * w;
                sumZ += z * w;
                sumW += w;
            }

            if (sumW == 0) return null;

            // Среднее и нормализация
            double mx = sumX / sumW;
            double my = sumY / sumW;
            double mz = sumZ / sumW;

            double hyp = Math.Sqrt(mx * mx + my * my);
            double latRes = Math.Atan2(mz, hyp);
            double lonRes = Math.Atan2(my, mx);

            // 5. Переводим обратно в Меркатор
            var (xRes, yRes) = SphericalMercator.FromLonLat(RadToDeg(lonRes), RadToDeg(latRes));
            return new MPoint(xRes, yRes);
        }

        private MPoint Project(double lon, double lat)
        {
            var (x, y) = SphericalMercator.FromLonLat(lon, lat);
            return new MPoint(x, y);
        }

        private LineString BuildRay(MPoint origin, double bearingDeg, double length = 200_000)
        {
            double angleRad = (90 - bearingDeg) * Math.PI / 180.0;
            var dx = length * Math.Cos(angleRad);
            var dy = length * Math.Sin(angleRad);

            var start = new Coordinate(origin.X, origin.Y);
            var end = new Coordinate(origin.X + dx, origin.Y + dy);
            return new LineString(new[] { start, end });
        }

        private static double DegToRad(double deg) => deg * Math.PI / 180.0;
        private static double RadToDeg(double rad) => rad * 180.0 / Math.PI;
    }
}
