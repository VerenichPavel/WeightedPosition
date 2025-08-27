using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using Mapsui;
using Mapsui.Layers;
using Mapsui.Nts;
using Mapsui.Providers;
using Mapsui.Styles;
using Mapsui.UI.Wpf;
using Mapsui.Tiling;
using Mapsui.Projections;
using NetTopologySuite.Geometries;
using System.Globalization;

// Алиасы для читаемости
using MapsuiBrush = Mapsui.Styles.Brush;
using NtsPoint = NetTopologySuite.Geometries.Point;
using Color = Mapsui.Styles.Color;
using Pen = Mapsui.Styles.Pen;

//  типы из DLL 
using WeightedPosition.Models;
using WeightedPosition.Services;
using ExCSS;

namespace TestApp.ViewModels
{   
    public partial class MainViewModel : ObservableObject
    {
        private readonly WeightedPositionCalculator _calculator = new();

        // Границы для авто-длины лучей (в метрах)
        private const double MinRayMeters = 1000.0;   
        private const double MaxRayMeters = 30000.0;  
        private const double RayViewportShare = 1; 

        [ObservableProperty]
        private string inputText =
@"A 53.9023 27.5615 60 0.9
B 53.9150 27.4500 120 0.8
C 53.9100 27.5800 300 0.7";

        [RelayCommand]
        private void Calculate(MapControl mapControl)
        {
            if (mapControl?.Map == null) return;

            var stations = ParseStations(InputText);

            // Пересобираем карту, сохраняем только подложку
            mapControl.Map.Layers.Clear();
            mapControl.Map.Layers.Add(OpenStreetMap.CreateTileLayer());

            if (stations.Count == 0) return;

            // Авто-длина лучей от текущего масштаба/размера
            var rayLengthMeters = GetDynamicRayLength(mapControl, stations);

            // 1) станции
            RenderStations(mapControl, stations);

            // 2) лучи (азимуты) + подписи с азимутом и длиной
            RenderBearings(mapControl, stations, MaxRayMeters);

            // 3) точки пересечений между лучами (локальная плоскостная аппроксимация)
            var intersections = CalculateIntersectionsPlanar(stations, rayLengthMeters * 4);
            RenderIntersections(mapControl, intersections);

            // 4) взвешенный центр из DLL
            var center = _calculator.CalculateWeightedCenter(stations);
            RenderCenter(mapControl, center);

            // 5) наводим камеру
            ZoomToContent(mapControl, stations, intersections, center);
        }

        [RelayCommand]
        private void Clear(MapControl mapControl)
        {
            InputText = string.Empty;

            if (mapControl?.Map == null) return;
            mapControl.Map.Layers.Clear();
            mapControl.Map.Layers.Add(OpenStreetMap.CreateTileLayer());
        }

        // --------------------------
        // Парсинг и нормализация
        // --------------------------

        private List<Station> ParseStations(string text)
        {
            var result = new List<Station>();
            if (string.IsNullOrWhiteSpace(text)) return result;

            var lines = text.Split(new[] { "\r\n", "\n" }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var raw in lines)
            {
                var parts = raw.Split(new[] { ' ', '\t', ';', ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length < 5) continue;

                if (!double.TryParse(parts[1], NumberStyles.Float, CultureInfo.InvariantCulture, out var lat)) continue;
                if (!double.TryParse(parts[2], NumberStyles.Float, CultureInfo.InvariantCulture, out var lon)) continue;
                if (!double.TryParse(parts[3], NumberStyles.Float, CultureInfo.InvariantCulture, out var brg)) continue;
                if (!double.TryParse(parts[4], NumberStyles.Float, CultureInfo.InvariantCulture, out var q)) continue;

                result.Add(new Station
                {
                    Name = parts[0],
                    LatitudeDeg = lat,
                    LongitudeDeg = lon,
                    BearingDeg = NormalizeDeg(brg),
                    Quality = Math.Clamp(q, 0.0, 10.0)
                });
            }
            return result;
        }

        private static double NormalizeDeg(double deg)
        {
            var a = deg % 360.0;
            return a < 0 ? a + 360.0 : a;
        }

        // --------------------------
        // Рендеринг
        // --------------------------

        private void RenderStations(MapControl mapControl, List<Station> stations)
        {
            var feats = stations.Select(s =>
            {
                var (x, y) = SphericalMercator.FromLonLat(s.LongitudeDeg, s.LatitudeDeg);
                var gf = new GeometryFeature { Geometry = new NtsPoint(x, y) };
                gf.Styles.Add(new SymbolStyle
                {
                    SymbolScale = 0.9,
                    Fill = new MapsuiBrush(Color.Red),
                    Outline = new Pen(Color.White, 1)
                });
                gf.Styles.Add(new LabelStyle
                {
                    Text = $"{s.Name} (Q={s.Quality:0.##})",
                    BackColor = new MapsuiBrush(new Color(255, 255, 255, 200)),
                    Halo = new Pen(Color.Black, 1),
                    Offset = new Offset(12, 0)
                });
                return gf;
            }).ToList();

            mapControl.Map.Layers.Add(new MemoryLayer
            {
                Name = "Stations",
                Features = feats
            });
        }

        private void RenderBearings(MapControl mapControl, List<Station> stations, double lengthMeters)
        {
            var lineFeats = new List<GeometryFeature>();
            var labelFeats = new List<GeometryFeature>();

            foreach (var s in stations)
            {
                // Начало (в метрах в проекции Меркатора)
                var start = SphericalMercator.FromLonLat(s.LongitudeDeg, s.LatitudeDeg);

                // Конец по азимуту на заданной длине (геодезический расчёт на сфере)
                var (lat2, lon2) = DestinationPoint(s.LatitudeDeg, s.LongitudeDeg, s.BearingDeg, lengthMeters);
                var end = SphericalMercator.FromLonLat(lon2, lat2);

                // Линия-луч
                var line = new LineString(new[]
                {
                    new Coordinate(start.x, start.y),
                    new Coordinate(end.x, end.y)
                });

                var lineGf = new GeometryFeature { Geometry = line };
                lineGf.Styles.Clear();
                lineGf.Styles.Add(new VectorStyle
                {
                    Line = new Pen(Color.Blue, 2)  // толщина 2, синий цвет
                });
                lineFeats.Add(lineGf);

                // Подпись вдоль луча (~70% длины)
                var labelX = start.x + 0.7 * (end.x - start.x);
                var labelY = start.y + 0.7 * (end.y - start.y);

                var labelGf = new GeometryFeature { Geometry = new NtsPoint(labelX, labelY) };
                labelGf.Styles.Add(new LabelStyle
                {
                    Text = $"{s.Name}: Az {s.BearingDeg:0.#}° • {lengthMeters / 1000.0:0.#} км",
                    BackColor = new MapsuiBrush(new Color(255, 255, 255, 220)),
                    Halo = new Pen(Color.Black, 1),
                    Offset = new Offset(0, 0)
                });
                labelFeats.Add(labelGf);
            }

            mapControl.Map.Layers.Add(new MemoryLayer
            {
                Name = "Bearings",
                Features = lineFeats
            });

            mapControl.Map.Layers.Add(new MemoryLayer
            {
                Name = "BearingLabels",
                Features = labelFeats
            });
        }

        private void RenderIntersections(MapControl mapControl, List<(double lat, double lon)> points)
        {
            if (points.Count == 0) return;

            var feats = points.Select(p =>
            {
                var (x, y) = SphericalMercator.FromLonLat(p.lon, p.lat);
                var gf = new GeometryFeature { Geometry = new NtsPoint(x, y) };
                gf.Styles.Add(new SymbolStyle
                {
                    SymbolScale = 0.8,
                    Fill = new MapsuiBrush(Color.Yellow),
                    Outline = new Pen(Color.Black, 1)
                });
                return gf;
            }).ToList();

            mapControl.Map.Layers.Add(new MemoryLayer
            {
                Name = "Intersections",
                Features = feats
            });
        }

        private void RenderCenter(MapControl mapControl, MPoint? center)
        {
            if (center is null) return;

            var feat = new GeometryFeature { Geometry = new NtsPoint(center.X, center.Y) };
            feat.Styles.Add(new SymbolStyle
            {
                SymbolScale = 1.2,
                Fill = new MapsuiBrush(Color.Black),
                Outline = new Pen(Color.White, 2)
            });
            feat.Styles.Add(new LabelStyle
            {
                Text = "Weighted",
                BackColor = new MapsuiBrush(new Color(255, 255, 255, 200)),
                Halo = new Pen(Color.Black, 1),
                Offset = new Offset(12, 0)
            });

            mapControl.Map.Layers.Add(new MemoryLayer
            {
                Name = "WeightedCenter",
                Features = new[] { feat }
            });
        }

        private void ZoomToContent(MapControl mapControl,
                                   List<Station> stations,
                                   List<(double lat, double lon)> intersections,
                                   MPoint? center)
        {
            // Собираем все точки в проекции карты, чтобы построить охватывающий BBox
            var points = new List<MPoint>();

            foreach (var s in stations)
            {
                var (x, y) = SphericalMercator.FromLonLat(s.LongitudeDeg, s.LatitudeDeg);
                points.Add(new MPoint(x, y));
            }

            foreach (var p in intersections)
            {
                var (x, y) = SphericalMercator.FromLonLat(p.lon, p.lat);
                points.Add(new MPoint(x, y));
            }

            if (center is not null) points.Add(center);

            if (points.Count == 0) return;

            var minX = points.Min(p => p.X);
            var minY = points.Min(p => p.Y);
            var maxX = points.Max(p => p.X);
            var maxY = points.Max(p => p.Y);

            // Немного полей
            var padding = 200.0;
            var bbox = new MRect(minX - padding, minY - padding, maxX + padding, maxY + padding);

            mapControl.Map.Navigator.ZoomToBox(bbox);
        }

        // --------------------------
        // Геометрия и геодезия
        // --------------------------

        // Длина луча как доля от текущей ширины видимой области карты (в метрах), с fallback'ами
        private double GetDynamicRayLength(MapControl mapControl, List<Station> stations)
        {
            // 1) Пытаемся взять текущую "разрешающую способность": метры на пиксель
            double? resolution = null;
            try
            {
                var resolutio = mapControl.Map.Navigator.Viewport.Resolution;
            }
            catch
            {
                // игнорируем — перейдём к запасным вариантам
            }

            if (resolution is double res && res > 0 && mapControl.ActualWidth > 0)
            {
                var viewportWidthMeters = res * mapControl.ActualWidth;
                var length = viewportWidthMeters * RayViewportShare;
                return Math.Clamp(length, MinRayMeters, MaxRayMeters);
            }

            // 2) Fallback: по охвату станций (ширина bbox в проекции Меркатора)
            if (stations.Count >= 2)
            {
                var merc = stations
                    .Select(s => SphericalMercator.FromLonLat(s.LongitudeDeg, s.LatitudeDeg))
                    .ToList();

                var minX = merc.Min(p => p.x);
                var maxX = merc.Max(p => p.x);
                var width = Math.Max(1.0, maxX - minX);
                var length = width * 0.25; // четверть охвата
                return Math.Clamp(length, MinRayMeters, MaxRayMeters);
            }

            // 3) Запасное значение
            return 5000.0;
        }

        // Геодезический "destination point" (WGS84-сфера)
        private static (double lat, double lon) DestinationPoint(double latDeg, double lonDeg, double bearingDeg, double distanceMeters)
        {
            const double R = 6378137.0; // радиус Земли (WGS84)
            var δ = distanceMeters / R;
            var θ = bearingDeg * Math.PI / 180.0;

            var φ1 = latDeg * Math.PI / 180.0;
            var λ1 = lonDeg * Math.PI / 180.0;

            var sinφ1 = Math.Sin(φ1);
            var cosφ1 = Math.Cos(φ1);
            var sinδ = Math.Sin(δ);
            var cosδ = Math.Cos(δ);
            var cosθ = Math.Cos(θ);
            var sinθ = Math.Sin(θ);

            var φ2 = Math.Asin(sinφ1 * cosδ + cosφ1 * sinδ * cosθ);
            var λ2 = λ1 + Math.Atan2(sinθ * sinδ * cosφ1, cosδ - sinφ1 * Math.Sin(φ2));

            var lat2 = φ2 * 180.0 / Math.PI;
            var lon2 = λ2 * 180.0 / Math.PI;
            return (lat2, lon2);
        }

        // Планарная аппроксимация пересечений двух лучей (для наглядности на масштабе ~5–50 км)
        private static List<(double lat, double lon)> CalculateIntersectionsPlanar(List<Station> stations, double maxDistanceMeters)
        {
            var result = new List<(double lat, double lon)>();
            if (stations.Count < 2) return result;

            // Переводим в метры (Web Mercator)
            var pts = stations.Select(s =>
            {
                var (x, y) = SphericalMercator.FromLonLat(s.LongitudeDeg, s.LatitudeDeg);
                // Единичный вектор направления по азимуту в локальной плоскости (ENU): X ~ East, Y ~ North
                var rad = s.BearingDeg * Math.PI / 180.0;
                var dir = new MPoint(Math.Sin(rad), Math.Cos(rad));
                return (s, x, y, dir);
            }).ToList();

            for (int i = 0; i < pts.Count; i++)
            {
                for (int j = i + 1; j < pts.Count; j++)
                {
                    var p = pts[i];
                    var q = pts[j];

                    // Решаем p0 + t*vp = q0 + s*vq, t>=0, s>=0
                    var p0x = p.x; var p0y = p.y;
                    var q0x = q.x; var q0y = q.y;

                    var vpx = p.dir.X; var vpy = p.dir.Y;
                    var vqx = q.dir.X; var vqy = q.dir.Y;

                    var denom = (vpx * vqy - vpy * vqx);
                    if (Math.Abs(denom) < 1e-9) continue; // почти параллельны

                    var dx = q0x - p0x;
                    var dy = q0y - p0y;

                    var t = (dx * vqy - dy * vqx) / denom;
                    var sParam = (dx * vpy - dy * vpx) / denom;

                    if (t < 0 || sParam < 0) continue;

                    var ix = p0x + t * vpx;
                    var iy = p0y + t * vpy;

                    // Ограничиваем по разумной дистанции от обеих станций
                    var dp = Math.Sqrt((ix - p0x) * (ix - p0x) + (iy - p0y) * (iy - p0y));
                    var dq2 = Math.Sqrt((ix - q0x) * (ix - q0x) + (iy - q0y) * (iy - q0y));
                    if (dp > maxDistanceMeters || dq2 > maxDistanceMeters) continue;

                    var (lon, lat) = SphericalMercator.ToLonLat(ix, iy);
                    result.Add((lat, lon));
                }
            }

            return result;
        }
    }
}

