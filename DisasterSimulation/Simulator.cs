using System;
using System.Collections.Generic;
using System.Drawing;
using System.Formats.Asn1;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json.Serialization;
using System.Threading.Tasks;
using DisasterSimulation;
using static System.Runtime.InteropServices.JavaScript.JSType;
using static DisasterSimulation.Simulator;

namespace DisasterSimulation
{
    public class Result
    {
        [JsonPropertyName("time")] public required double Time { get; set; }
        [JsonPropertyName("particlePositions")] public required Vector3?[] ParticlePositions { get; set; }
    }

    public class Setting
    {
        [JsonPropertyName("influenceRadius")]
        public required double InfluenceRadius { get; set; } //影響半径

        [JsonPropertyName("particleMass")]
        public required double ParticleMass { get; set; } //粒子の質量


        [JsonPropertyName("gravitationalAcceleration")]
        public required  Vector3 GravitationalAcceleration { get; set; }  // 重力加速度

        [JsonPropertyName("pressureStiffness")]
        public required double PressureStiffness { get; set; } //圧力係数

        [JsonPropertyName("restDensity")]
        public required double RestDensity { get; set; } //静止密度

        [JsonPropertyName("viscosity")]
        public required double Viscosity { get; set; }  // 粘性係数

        [JsonPropertyName("dampingCoefficient")]
        public required double DampingCoefficient { get; set; }  // ダンパ係数

        [JsonPropertyName("springConstant")]
        public required double SpringConstant { get; set; }  // ばね係数
    }
    
    internal class FaceData
    {
        [JsonPropertyName("index")] public int Index { get; set; }
        [JsonPropertyName("triangle")] public required Vector3[] Triangle { get; set; }
        [JsonPropertyName("centerOfGravity")] public required Vector3 CenterOfGravity { get; set; }
        [JsonPropertyName("normalVector")] public required Vector3 NormalVector { get; set; }
        [JsonPropertyName("insideJudge")] public required List<InsideJudge> InsideJudges { get; set; }
        public class InsideJudge
        {
            [JsonPropertyName("point")] public required Vector3 Point { get; set; }
            [JsonPropertyName("normalVector")] public required Vector3 NormalVector { get; set; }
        }
    }

    internal class FaceDataInfo
    {
        public double maxX = double.NegativeInfinity;
        public double minX = double.PositiveInfinity;
        public double maxY = double.NegativeInfinity;
        public double minY = double.PositiveInfinity;
        public double maxZ = double.NegativeInfinity;
        public double minZ = double.PositiveInfinity;
    }

    internal class Simulator(FaceData[] data)
    {
        public class Particle
        {
            public Vector3 position = new();
            public Vector3 velocity = new();
            public Vector3 force = new();
            public Vector3 acceleration = new();
            public double density;
            public double pressure;
            public uint id;
            public List<uint> affectingParticleIds = new();
            public List<double> DistancesBetweenAffectingParticle = new();
            public List<Vector3> vectorsBetweenAffectingParticle = new();
            public List<int> affectingParticleIndex = new();
        }
        public class ParticlePositionInSingleDirectionWithId
        {
            public double pos;
            public uint id;
        }
        class Term
        {
            public Vector3 pressureTerm = new();
            public Vector3 viscosityTerm = new();
            public Vector3 coliderTerm = new();
        }

        internal List<Result> result = [];

        readonly FaceData[] data = data;
        FaceDataInfo[] dataInfo = new FaceDataInfo[data.Length];

        readonly Dictionary<int, Term> terms = [];

        readonly static double tsunamiSpeed = 10;

        readonly static double numberOfNearParticle = 12;

        readonly static double particleDistance = 2/*0.01*/;
        readonly static double volumePerParticle = Math.Pow(particleDistance, 3d);
        readonly static double particleDiameter = Math.Pow(volumePerParticle, 1d / 3d);

        readonly static Vector3 g = new(0, -9.8, 0);  // 重力加速度
        readonly static double pressureStiffness = 0.1/*200*/; //圧力係数
        readonly static double restDensity = 1000; //静止密度
        readonly static double viscosity = 0.000001;  // 粘性係数
        readonly static double dampingCoefficient = /*0.256*//*256*/0.03;  // ダンパ係数
        readonly static double springConstant = /*10*//*10000*//*7*/15;  // ばね係数

        readonly static double h = Math.Pow((3 * numberOfNearParticle) / (4 * Math.PI), 1d / 3d) * particleDiameter/*100*//*0.3*//* 0.012 */; //影響半径
        readonly static double particleMass = restDensity * volumePerParticle/*9000*//*0.0002*/; //粒子の質量

        //readonly double densityCoef = particleMass * 315 / (64 * Math.PI * Math.Pow(h, 9)); //密度計算で使う
        readonly double densityCoef = particleMass * 15 / (Math.PI * Math.Pow(h, 6)); //密度計算で使う
        readonly double pressureCoef = particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //圧力項計算で使う
        readonly double viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //粘性項計算で使う

        readonly List<Particle> _particles = [];
        public uint lastUsedId = 0;
        public List<ParticlePositionInSingleDirectionWithId> _particleXwithId = new();
        public List<ParticlePositionInSingleDirectionWithId> _particleYwithId = new();
        public List<ParticlePositionInSingleDirectionWithId> _particleZwithId = new();
        readonly static double deltaTime = 0.1/*0.03*/;

        Task CalcAffectingParticles(double h, List<Particle> particles, List<ParticlePositionInSingleDirectionWithId> particleXwithId, List<ParticlePositionInSingleDirectionWithId> particleYwithId, List<ParticlePositionInSingleDirectionWithId> particleZwithId)
        {
            List<Task> tasks = new List<Task>();
            for (int index = 0; index < particles.Count; index++)
            {
                int i = index;
                Task task = Task.Run(() =>
                {
                    uint selfId;
                    double mainNumInSearching;
                    int selectedIndex = 0;
                    int UpperLimitIndex;
                    int LowerLimitIndex;
                    List<ParticlePositionInSingleDirectionWithId> selectedYParticleDatas = new();
                    List<ParticlePositionInSingleDirectionWithId> selectedZParticleDatas = new();
                    selfId = particles[i].id;
                    List<uint> selectedIdList = new List<uint>();

                    //影響境界が粒子群の内部にあることを確認する
                    if (particles[i].position.X - h > particleXwithId[particleXwithId.Count - 1].pos)
                    {

                    }
                    else if (particles[i].position.X + h < particleXwithId[0].pos)
                    {

                    }
                    else if (particles[i].position.Y - h > particleYwithId[particleYwithId.Count - 1].pos)
                    {

                    }
                    else if (particles[i].position.Y + h < particleYwithId[0].pos)
                    {

                    }
                    else if (particles[i].position.Z - h > particleZwithId[particleZwithId.Count - 1].pos)
                    {

                    }
                    else if (particles[i].position.Z + h < particleZwithId[0].pos)
                    {

                    }
                    else
                    {
                        //X軸での下限を求める
                        mainNumInSearching = particles[i].position.X - h;
                        UpperLimitIndex = particleXwithId.Count - 1;
                        LowerLimitIndex = 0;
                        int smallestParticleInAffect_Index;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleXwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleXwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleXwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //X軸での上限を求める
                        mainNumInSearching = particles[i].position.X + h;
                        UpperLimitIndex = particleXwithId.Count - 1;
                        LowerLimitIndex = 0;
                        int largestParticleInAffect_Index;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleXwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleXwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleXwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        for (int j = smallestParticleInAffect_Index; j <= largestParticleInAffect_Index; j++)
                        {
                            selectedIdList.Add(particleXwithId[j].id);
                        }
                        for (int j = 0; j < selectedIdList.Count; j++)
                        {
                            if (selectedIdList[j] == selfId)
                            {
                                selectedIdList.RemoveAt(j);
                            }
                        }


                        //Y軸での下限を求める
                        mainNumInSearching = particles[i].position.Y - h;
                        UpperLimitIndex = particleYwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleYwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleYwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleYwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //Y軸での上限を求める
                        mainNumInSearching = particles[i].position.Y + h;
                        UpperLimitIndex = particleYwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleYwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleYwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleYwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        int Loop_containListRemove_Adder = 1;
                        bool alsoSelected = false;
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {
                            alsoSelected = false;
                            Loop_containListRemove_Adder = 1;
                            for (int k = smallestParticleInAffect_Index; k <= largestParticleInAffect_Index; k++)
                            {
                                if (selectedIdList[j] == particleYwithId[k].id)
                                {
                                    alsoSelected = true;
                                    break;
                                }
                            }
                            if (!alsoSelected)
                            {
                                selectedIdList.RemoveAt(j);
                                Loop_containListRemove_Adder = 0;
                            }
                        }

                        //Z軸での下限を求める
                        mainNumInSearching = particles[i].position.Z - h;
                        UpperLimitIndex = particleZwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleZwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleZwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleZwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //Z軸での上限を求める
                        mainNumInSearching = particles[i].position.Z + h;
                        UpperLimitIndex = particleZwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleZwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleZwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleZwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {
                            alsoSelected = false;
                            Loop_containListRemove_Adder = 1;
                            for (int k = smallestParticleInAffect_Index; k <= largestParticleInAffect_Index; k++)
                            {
                                if (selectedIdList[j] == particleZwithId[k].id)
                                {
                                    alsoSelected = true;
                                    break;
                                }
                            }
                            if (!alsoSelected)
                            {
                                selectedIdList.RemoveAt(j);
                                Loop_containListRemove_Adder = 0;
                            }
                        }
                        List<double> distancesBetweenSelectedParticle = new();
                        List<int> selectedParticleIndex = new();
                        List<Vector3> vectorsBetweenAffectingParticle = new();
                        double h2 = h * h;
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {

                            Loop_containListRemove_Adder = 1;
                            //particleはAddParticleを通して増加する場合、常にidが小さい順に並んでいて、必ず現在注目している粒子が存在することが保証されているので、二分探索を用いる
                            selectedIndex = 0;
                            UpperLimitIndex = particles.Count;
                            LowerLimitIndex = 0;
                            while (true)
                            {
                                selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                                if (selectedIdList[j] == particles[selectedIndex].id)
                                {
                                    break;
                                }
                                else if (selectedIdList[j] < particles[selectedIndex].id)
                                {
                                    UpperLimitIndex = selectedIndex;
                                }
                                else
                                {
                                    LowerLimitIndex = selectedIndex;
                                }
                            }


                            Vector3 diff = Vector3Utility.SubVector3(particles[selectedIndex].position, particles[i].position); //粒子距離
                            //これってベクトルの向きが注目粒子から影響半径内の粒子へ　になってるけど大丈夫なのか
                            double r2 = Vector3Utility.DotVector3(diff, diff); //粒子距離の２乗
                            if (r2 >= h2 || r2 == 0)
                            {
                                selectedIdList.Remove(selectedIdList[j]);
                                Loop_containListRemove_Adder = 0;
                            }
                            else
                            {
                                selectedParticleIndex.Add(selectedIndex);
                                distancesBetweenSelectedParticle.Add(Math.Sqrt(r2));
                                vectorsBetweenAffectingParticle.Add(diff);
                            }
                        }
                        particles[i].affectingParticleIds = selectedIdList;
                        particles[i].affectingParticleIndex = selectedParticleIndex;
                        particles[i].DistancesBetweenAffectingParticle = distancesBetweenSelectedParticle;
                        particles[i].vectorsBetweenAffectingParticle = vectorsBetweenAffectingParticle;
                    }
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);

        }

        Task CalcDensity(List<Particle> particles)
        {
            List<Task> tasks = [];
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {
                    Particle nowParticle = particles[index]; //今回計算する粒子
                    double sum = 0; //足し合わせる変数
                    for (int j = 0; j < nowParticle.affectingParticleIds.Count; j++)
                    {
                        double r2 = nowParticle.DistancesBetweenAffectingParticle[j] * nowParticle.DistancesBetweenAffectingParticle[j];
                        double c = h2 - r2;
                        sum += Math.Pow(c, 3);
                    }
                    
                    nowParticle.density = sum * densityCoef; //密度が求まった
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        static Task CalcPressure(List<Particle> particles)
        {
            List<Task> tasks = [];
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の圧力を計算
                int index = i;
                Task task = Task.Run(() =>
                {
                    double pressure = pressureStiffness * (particles[index].density - restDensity);
                    if (pressure >= 0)
                    {
                        particles[index].pressure = pressure;
                    }
                    else
                    {
                        particles[index].pressure = 0;
                    }
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Task CalcPressureTerm(List<Particle> particles)
        {
            List<Task> tasks = [];
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {
                    Particle nowParticle = particles[index]; //今回計算する粒子
                    //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                    Vector3 sum = new(); //足し合わせる変数
                    for (int j = 0; j < nowParticle.affectingParticleIds.Count; j++)
                    {
                        double r = nowParticle.DistancesBetweenAffectingParticle[j];
                        double c = h - r;
                        double n = ((particles[nowParticle.affectingParticleIndex[j]].pressure /*-*/
                        +nowParticle.pressure) / (2 * particles[nowParticle.affectingParticleIndex[j]].density)) * Math.Pow(c, 2) / r;
                        sum = Vector3Utility.AddVector3(sum, Vector3Utility.MultiplyScalarVector3(nowParticle.vectorsBetweenAffectingParticle[j], n));
                    }
                    terms[index].pressureTerm = Vector3Utility.MultiplyScalarVector3(sum, (-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
                    //Console.WriteLine($"{terms[index].pressureTerm.X},{terms[index].pressureTerm.Y},{terms[index].pressureTerm.Z}");
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Task CalcViscosityTerm(List<Particle> particles)
        {
            List<Task> tasks = [];
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {
                    Particle nowParticle = particles[index]; //今回計算する粒子
                    //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                    Vector3 sum = new(); //足し合わせる変数
                    for (int j = 0; j < nowParticle.affectingParticleIds.Count; j++)
                    {
                        double r = nowParticle.DistancesBetweenAffectingParticle[j];
                        double c = h - r;
                        double n = c / particles[nowParticle.affectingParticleIndex[j]].density;
                        sum = Vector3Utility.AddVector3(sum, Vector3Utility.MultiplyScalarVector3(Vector3Utility.SubVector3( particles[nowParticle.affectingParticleIndex[j]].velocity, nowParticle.velocity), n));
                    }

                    terms[index].viscosityTerm = Vector3Utility.MultiplyScalarVector3(sum, viscosityCoef);  // 粘性項が求まった
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        static bool Is_inside(Vector3 position, FaceData point)
        {
            for (int i = 0; i < point.InsideJudges.Count; i++)
            {
                double dot = Vector3Utility.DotVector3(point.InsideJudges[i].NormalVector, Vector3Utility.SubVector3(position, point.InsideJudges[i].Point));
                if (dot < 0) return false;
            }
            return true;
            /*const dot = dotVector3(point.normalVector, subVector3(position, point.centerOfGravity));
            return dot <= 0;*/
        }
        Task CalcColiderTerm(List<Particle> particles, List<ParticlePositionInSingleDirectionWithId> particleXwithId, List<ParticlePositionInSingleDirectionWithId> particleYwithId, List<ParticlePositionInSingleDirectionWithId> particleZwithId)
        {
            List<Task> tasks = [];
            //前回の演算の引き継ぎを防ぐために、termのColliderTermを初期化
            for (int i = 0; i < terms.Count; i++)
            {
                terms[i].coliderTerm = new Vector3();
            }

            for(int i = 0; i < data.Length; i++)
            {
                int index = i;
                Task task = Task.Run(() =>
                {
                    double mainNumInSearching;
                    int selectedIndex = 0;
                    int UpperLimitIndex;
                    int LowerLimitIndex;
                    List<ParticlePositionInSingleDirectionWithId> selectedYParticleDatas = new();
                    List<ParticlePositionInSingleDirectionWithId> selectedZParticleDatas = new();
                    List<uint> selectedIdList = new List<uint>();

                    //影響境界が粒子群の内部にあることを確認する
                    if (dataInfo[index].minX > particleXwithId[particleXwithId.Count - 1].pos)
                    {

                    }
                    else if (dataInfo[index].maxX < particleXwithId[0].pos)
                    {

                    }
                    else if (dataInfo[index].minY > particleYwithId[particleYwithId.Count - 1].pos)
                    {

                    }
                    else if (dataInfo[index].maxY < particleYwithId[0].pos)
                    {

                    }
                    else if (dataInfo[index].minZ > particleZwithId[particleZwithId.Count - 1].pos)
                    {

                    }
                    else if (dataInfo[index].maxZ < particleZwithId[0].pos)
                    {

                    }
                    else
                    {

                        //X軸での下限を求める
                        mainNumInSearching = dataInfo[index].minX;
                        UpperLimitIndex = particleXwithId.Count - 1;
                        LowerLimitIndex = 0;
                        int smallestParticleInAffect_Index;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleXwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleXwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleXwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //X軸での上限を求める
                        mainNumInSearching = dataInfo[index].maxX;
                        UpperLimitIndex = particleXwithId.Count - 1;
                        LowerLimitIndex = 0;
                        int largestParticleInAffect_Index;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleXwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleXwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleXwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        for (int j = smallestParticleInAffect_Index; j <= largestParticleInAffect_Index; j++)
                        {
                            selectedIdList.Add(particleXwithId[j].id);
                        }
                        //Y軸での下限を求める
                        mainNumInSearching = dataInfo[index].minY;
                        UpperLimitIndex = particleYwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleYwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleYwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleYwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //Y軸での上限を求める
                        mainNumInSearching = dataInfo[index].maxY;
                        UpperLimitIndex = particleYwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleYwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleYwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleYwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        int Loop_containListRemove_Adder = 1;
                        bool alsoSelected = false;
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {
                            alsoSelected = false;
                            Loop_containListRemove_Adder = 1;
                            for (int k = smallestParticleInAffect_Index; k <= largestParticleInAffect_Index; k++)
                            {
                                if (selectedIdList[j] == particleYwithId[k].id)
                                {
                                    alsoSelected = true;
                                    break;
                                }
                            }
                            if (!alsoSelected)
                            {
                                selectedIdList.RemoveAt(j);
                                Loop_containListRemove_Adder = 0;
                            }
                        }

                        //Z軸での下限を求める
                        mainNumInSearching = dataInfo[index].minZ;
                        UpperLimitIndex = particleZwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleZwithId[LowerLimitIndex].pos >= mainNumInSearching)
                                {
                                    smallestParticleInAffect_Index = LowerLimitIndex;
                                }
                                else
                                {
                                    smallestParticleInAffect_Index = UpperLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleZwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                smallestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleZwithId[selectedIndex].pos < mainNumInSearching)
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                            else
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                        }
                        //Z軸での上限を求める
                        mainNumInSearching = dataInfo[index].maxZ;
                        UpperLimitIndex = particleZwithId.Count - 1;
                        LowerLimitIndex = 0;
                        while (true)
                        {
                            if (UpperLimitIndex - LowerLimitIndex <= 1)
                            {
                                //探索終了
                                if (particleZwithId[UpperLimitIndex].pos <= mainNumInSearching)
                                {
                                    largestParticleInAffect_Index = UpperLimitIndex;
                                }
                                else
                                {
                                    largestParticleInAffect_Index = LowerLimitIndex;
                                }
                                break;
                            }
                            selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                            if (particleZwithId[selectedIndex].pos == mainNumInSearching)
                            {
                                largestParticleInAffect_Index = selectedIndex;
                                break;
                            }
                            else if (particleZwithId[selectedIndex].pos > mainNumInSearching)
                            {
                                UpperLimitIndex = selectedIndex;
                            }
                            else
                            {
                                LowerLimitIndex = selectedIndex;
                            }
                        }
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {
                            alsoSelected = false;
                            Loop_containListRemove_Adder = 1;
                            for (int k = smallestParticleInAffect_Index; k <= largestParticleInAffect_Index; k++)
                            {
                                if (selectedIdList[j] == particleZwithId[k].id)
                                {
                                    alsoSelected = true;
                                    break;
                                }
                            }
                            if (!alsoSelected)
                            {
                                selectedIdList.RemoveAt(j);
                                Loop_containListRemove_Adder = 0;
                            }
                        }
                        //idからindexを特定し、内在判定を行い、内在していれば反発力を計算する
                        for (int j = 0; j < selectedIdList.Count; j += Loop_containListRemove_Adder)
                        {
                            Loop_containListRemove_Adder = 1;
                            selectedIndex = 0;
                            UpperLimitIndex = particles.Count;
                            LowerLimitIndex = 0;
                            while (true)
                            {
                                selectedIndex = (UpperLimitIndex + LowerLimitIndex) / 2;
                                if (selectedIdList[j] == particles[selectedIndex].id)
                                {
                                    break;
                                }
                                else if (selectedIdList[j] < particles[selectedIndex].id)
                                {
                                    UpperLimitIndex = selectedIndex;
                                }
                                else
                                {
                                    LowerLimitIndex = selectedIndex;
                                }
                            }
                            Vector3 term = new();
                            //内在判定を行う
                            if (Is_inside(particles[selectedIndex].position, data[index]))
                            {
                                double distance = Math.Abs(Vector3Utility.DotVector3(Vector3Utility.SubVector3(particles[selectedIndex].position, data[index].CenterOfGravity), data[index].NormalVector));
                                Vector3 nowTerm = Vector3Utility.MultiplyScalarVector3(data[index].NormalVector, springConstant * distance + dampingCoefficient * Vector3Utility.DotVector3(particles[selectedIndex].velocity, data[index].NormalVector));
                                term = Vector3Utility.AddVector3(term, nowTerm);
                            }
                            terms[selectedIndex].coliderTerm = term;
                        }
                    }
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Vector3?[] Tick()
        {
            for (int i = 0; i < _particles.Count; i++)
            {
                if (!terms.ContainsKey(i))
                {
                    terms.Add(i, new Term());
                }
            }
            //粒子ごとに影響を与える粒子の配列を求める
            CalcAffectingParticles(h,_particles, _particleXwithId, _particleYwithId, _particleZwithId).Wait();
            CalcDensity(_particles).Wait();
            CalcPressure(_particles).Wait();
            CalcPressureTerm(_particles).Wait();
            CalcViscosityTerm(_particles).Wait();
            CalcColiderTerm(_particles,_particleXwithId,_particleYwithId,_particleZwithId).Wait();

            _particleXwithId.Clear();
            _particleYwithId.Clear();
            _particleZwithId.Clear();

            List<Vector3?> tickResult = [];
            for (int i = 0; i < _particles.Count; i++)
            {
                Particle nowParticle = _particles[i];
                //if (nowParticle.is_wall) continue;
                //粒子位置の計算
                Vector3 a;

                a = Vector3Utility.AddVector3(Vector3Utility.AddVector3(Vector3Utility.AddVector3(terms[i].pressureTerm, terms[i].viscosityTerm), terms[i].coliderTerm), g);
                Vector3 v = Vector3Utility.AddVector3(nowParticle.velocity, Vector3Utility.MultiplyScalarVector3(Vector3Utility.AddVector3(nowParticle.acceleration, a), 0.5 * deltaTime));
                Vector3 deltaPosition = Vector3Utility.AddVector3(Vector3Utility.MultiplyScalarVector3(v, deltaTime), Vector3Utility.MultiplyScalarVector3(a, 0.5 * deltaTime * deltaTime));
                nowParticle.acceleration = a;
                nowParticle.velocity = v;
                nowParticle.position = Vector3Utility.AddVector3(nowParticle.position, deltaPosition);
                
                tickResult.Add(nowParticle.position);
            }
            _particles.RemoveAll(particle => particle.position.Y <= 0);
            return [.. tickResult];
        }

        static void AddParticles(List<Particle> particles, ref uint lastUsedId)
        {
            Vector3 tsunamiVelocity = Vector3Utility.MultiplyScalarVector3(Vector3Utility.NormalizeVector3(new Vector3(0.72342719346605268235776097708833, 0, 1)), tsunamiSpeed);
            for (double z = -1288; z <= 1155; z += particleDistance)
            {
                for (double y = 38; y <= 43; y += particleDistance)
                {
                    Particle particle = new()
                    {
                        position = new Vector3(0.72342719346605268235776097708833 * z - 582.38286403266973658821119511456, y, z),
                        velocity = tsunamiVelocity,
                        id = lastUsedId
                    };
                    lastUsedId = lastUsedId + 1;
                    particles.Add(particle);
                }
            }
        }
        static void SortParticlesOnStart(List<Particle> particles, List<ParticlePositionInSingleDirectionWithId> PPsISDWId, string direction)
        {
            for (int i = 0; i < particles.Count; i++)
            {
                double particlePos;
                if (direction == "x")
                {
                    particlePos = particles[i].position.X;
                }
                else if (direction == "y")
                {
                    particlePos = particles[i].position.Y;
                }
                else
                {
                    particlePos = particles[i].position.Z;
                }
                ParticlePositionInSingleDirectionWithId PPISDWId = new()
                {
                    pos = particlePos,
                    id = particles[i].id
                };
                PPsISDWId.Add(PPISDWId);
            }
            int unSearchedRange = particles.Count;
            double largestPos;
            int largestPos_index;
            uint largestPos_id;
            for (int i = 0; i < particles.Count; i++)
            {
                largestPos = PPsISDWId[0].pos;
                largestPos_index = 0;
                largestPos_id = PPsISDWId[0].id;
                for (int j = 0; j < unSearchedRange; j++)
                {
                    
                    if (PPsISDWId[j].pos >= largestPos)
                    {
                        largestPos = PPsISDWId[j].pos;
                        largestPos_index = j;
                        largestPos_id = PPsISDWId[j].id;
                    }
                }
                PPsISDWId[largestPos_index].pos = PPsISDWId[unSearchedRange - 1].pos;
                PPsISDWId[largestPos_index].id = PPsISDWId[unSearchedRange - 1].id;
                PPsISDWId[unSearchedRange - 1].pos = largestPos;
                PPsISDWId[unSearchedRange - 1].id = largestPos_id;
                unSearchedRange--;
            }

        }
        public void Start(double simulateSeconds)
        {
            // 初めに実行する処理
            //makeWall(_particles);

            Console.WriteLine("シミュレーション開始");

            double height = 5;  // 最大食い込み量（食い込んでいる面を特定するために用いる）←createFaces.jsと同じ値を用いること
            for (int i = 0; i < data.Length; i++)
            {
                dataInfo[i] = new();
                FaceData point = data[i];
                Vector3 heightVector = Vector3Utility.MultiplyScalarVector3(point.NormalVector, -height);
                for (int j = 0; j < point.Triangle.Length; j++)
                {
                    Vector3 point1 = point.Triangle[j];
                    Vector3 point2 = Vector3Utility.AddVector3(point.Triangle[j], heightVector);
                    if (point1.X < dataInfo[i].minX)
                    {
                        dataInfo[i].minX = point1.X;
                    }
                    if (point1.X > dataInfo[i].maxX)
                    {
                        dataInfo[i].maxX = point1.X;
                    }
                    if (point1.Y < dataInfo[i].minY)
                    {
                        dataInfo[i].minY = point1.Y;
                    }
                    if (point1.Y > dataInfo[i].maxY)
                    {
                        dataInfo[i].maxY = point1.Y;
                    }
                    if (point1.Z < dataInfo[i].minZ)
                    {
                        dataInfo[i].minZ = point1.Z;
                    }
                    if (point1.Z > dataInfo[i].maxZ)
                    {
                        dataInfo[i].maxZ = point1.Z;
                    }
                    if (point2.X < dataInfo[i].minX)
                    {
                        dataInfo[i].minX = point2.X;
                    }
                    if (point2.X > dataInfo[i].maxX)
                    {
                        dataInfo[i].maxX = point2.X;
                    }
                    if (point2.Y < dataInfo[i].minY)
                    {
                        dataInfo[i].minY = point2.Y;
                    }
                    if (point2.Y > dataInfo[i].maxY)
                    {
                        dataInfo[i].maxY = point2.Y;
                    }
                    if (point2.Z < dataInfo[i].minZ)
                    {
                        dataInfo[i].minZ = point2.Z;
                    }
                    if (point2.Z > dataInfo[i].maxZ)
                    {
                        dataInfo[i].maxZ = point2.Z;
                    }
                }
            }
            
            for (int i = 0; i < simulateSeconds / deltaTime; i++)
            {
                if (i % 3  == 0)
                {
                    AddParticles(_particles, ref lastUsedId);
                }
                SortParticlesOnStart(_particles, _particleXwithId, "x");
                SortParticlesOnStart(_particles, _particleYwithId, "y");
                SortParticlesOnStart(_particles, _particleZwithId, "z");

                
                double time = (i + 1) * deltaTime;
                string reportHeader = time + "秒目/" + simulateSeconds + "秒 ";
                Console.WriteLine(reportHeader);
                Vector3?[] particlePositions = Tick();
                Result tickResult = new()
                {
                    Time = time,
                    ParticlePositions = particlePositions
                };
                result.Add(tickResult);
            }
            
        }
    }
}
