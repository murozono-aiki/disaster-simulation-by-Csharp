using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.Json.Serialization;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace DisasterSimulation
{
    internal class Vector3(double x = 0, double y = 0, double z = 0)
    {
        [JsonPropertyName("x")] public double X { get; set; } = x;
        [JsonPropertyName("y")] public double Y { get; set; } = y;
        [JsonPropertyName("z")] public double Z { get; set; } = z;
    }
    internal static class Vector3Utility
    {
        public static Vector3 AddVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }

        public static Vector3 SubVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Vector3 MultiplyScalarVector3(Vector3 a, double n)
        {
            return new Vector3(n * a.X, n * a.Y, n * a.Z);
        }

        public static double DotVector3(Vector3 a, Vector3 b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        }

        public static Vector3 CrossVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
        }
    }
}
