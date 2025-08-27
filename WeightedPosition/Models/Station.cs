using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WeightedPosition.Models
{
    public class Station
    {
        public string Name { get; set; } = string.Empty;
        public double LatitudeDeg { get; set; }
        public double LongitudeDeg { get; set; }
        public double BearingDeg { get; set; } // навигационный азимут
        public double Quality { get; set; } = 1.0; // вес пеленга
    }
}
