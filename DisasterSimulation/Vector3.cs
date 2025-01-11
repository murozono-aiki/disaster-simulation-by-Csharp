using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace DisasterSimulation
{
    internal class Vector3(double x = 0, double y = 0, double z = 0)
    {
        public double x = x;
        public double y = y;
        public double z = z;
    }
    internal static class Vector3Utility
    {
        public static Vector3 AddVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static Vector3 SubVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static Vector3 MultiplyScalarVector3(Vector3 a, double n)
        {
            return new Vector3(n * a.x, n * a.y, n * a.z);
        }

        public static double DotVector3(Vector3 a, Vector3 b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }

        public static Vector3 CrossVector3(Vector3 a, Vector3 b)
        {
            return new Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }
    }
}
