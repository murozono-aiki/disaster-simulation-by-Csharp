using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace DisasterSimulation
{
    internal class Simulator
    {
    }

    internal class PointData
    {
        public int index;
        public Vector3[] triangle = new Vector3[3];
        public Vector3 centerOfGravity;
        public Vector3 normalVector;
        public InsideJudge[] insideJudge;
        public class InsideJudge
        {
            public Vector3 point;
            public Vector3 normalVector;
        }
    }
}
