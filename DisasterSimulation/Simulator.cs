using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics;
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
        [JsonPropertyName("h")]
        public required double H { get; set; } //影響半径

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

        [JsonPropertyName("attenuationCoefficient")]
        public required double AttenuationCoefficient { get; set; }  // ダンパ係数

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
        }
        class Term
        {
            public Vector3 pressureTerm = new();
            public Vector3 viscosityTerm = new();
            public Vector3 coliderTerm = new();
        }

        internal List<Result> result = new();

        readonly FaceData[] data = data;

        Dictionary<int, Term> terms = new();

        readonly static double h = 0.3/* 0.012 */; //影響半径
        readonly static double particleMass = 0.0002; //粒子の質量

        readonly static Vector3 g = new(0, -9.8, 0);  // 重力加速度
        readonly static double pressureStiffness = 200; //圧力係数
        readonly static double restDensity = 1000; //静止密度
        readonly static double viscosity = 0.000001;  // 粘性係数
        readonly static double attenuationCoefficient = -5;  // ダンパ係数（javascript側と同じ値を用いること）
        readonly static double springConstant = -5;  // ばね係数

        readonly double densityCoef = particleMass * 315 / (64 * Math.PI * Math.Pow(h, 9)); //密度計算で使う

        readonly double pressureCoef = particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //圧力項計算で使う
        readonly double viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //粘性項計算で使う

        List<Particle> _particles = new();
        readonly static double deltaTime = 0.1;

        Task CalcDensity(List<Particle> particles)
        {
            List<Task> tasks = new();
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {

                    //reportProgress(reportHeader + "密度を計算中..." + i + "/" + particles.Count);

                    Particle nowParticle = particles[index]; //今回計算する粒子
                    double sum = 0; //足し合わせる変数
                    for (int j = 0; j < particles.Count; j++)
                    { //他の粒子全てについて
                        if (index == j) { continue; } //自分自身だったらスキップ
                        Particle nearParticle = particles[j];

                        Vector3 diff = Vector3Utility.SubVector3(nearParticle.position, nowParticle.position); //粒子距離
                        double r2 = Vector3Utility.DotVector3(diff, diff); //粒子距離の２乗

                        //粒子距離がhより小さい場合だけ計算する
                        if (r2 < h2)
                        {
                            double c = h2 - r2;
                            sum += Math.Pow(c, 3); //(h2-r2)の３乗
                        }
                    }

                    nowParticle.density = sum * densityCoef; //密度が求まった
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Task CalcPressure(List<Particle> particles)
        {
            List<Task> tasks = new();
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の圧力を計算
                int index = i;
                Task task = Task.Run(() =>
                {

                    //reportProgress(reportHeader + "圧力を計算中..." + i + "/" + particles.Count);

                    particles[index].pressure = pressureStiffness * (particles[index].density - restDensity);
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Task CalcPressureTerm(List<Particle> particles)
        {
            List<Task> tasks = new();
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {

                    //reportProgress(reportHeader + "圧力項を計算中..." + i + "/" + particles.length);

                    Particle nowParticle = particles[index]; //今回計算する粒子
                    //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                    Vector3 sum = new(); //足し合わせる変数
                    for (int j = 0; j < particles.Count; j++)
                    { //他の粒子全てについて
                        if (index == j) { continue; } //自分自身だったらスキップ
                        Particle nearParticle = particles[j];

                        Vector3 diff = Vector3Utility.SubVector3(nearParticle.position, nowParticle.position); //粒子距離
                        double r2 = Vector3Utility.DotVector3(diff, diff); //粒子距離の２乗

                        //粒子距離がhより小さい場合だけ計算する
                        if (r2 < h2)
                        {
                            double r = Math.Sqrt(r2); //粒子距離
                            double c = h - r;
                            double n = ((nearParticle.pressure /*-*/+ nowParticle.pressure) / (2 * nearParticle.density)) * Math.Pow(c, 2) / r;
                            sum = Vector3Utility.AddVector3(sum, Vector3Utility.MultiplyScalarVector3(diff, n));
                        }
                    }

                    terms[index].pressureTerm = Vector3Utility.MultiplyScalarVector3(sum, (-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
                });
                tasks.Add(task);
            }
            return Task.WhenAll(tasks);
        }

        Task CalcViscosityTerm(List<Particle> particles)
        {
            List<Task> tasks = new();
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算
                int index = i;
                Task task = Task.Run(() =>
                {

                    //reportProgress(reportHeader + "粘性項を計算中..." + i + "/" + particles.length);

                    Particle nowParticle = particles[index]; //今回計算する粒子
                    //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                    Vector3 sum = new(); //足し合わせる変数
                    for (int j = 0; j < particles.Count; j++)
                    { //他の粒子全てについて
                        if (index == j) { continue; } //自分自身だったらスキップ
                        Particle nearParticle = particles[j];

                        Vector3 diff = Vector3Utility.SubVector3(nearParticle.position, nowParticle.position); //粒子距離
                        double r2 = Vector3Utility.DotVector3(diff, diff); //粒子距離の２乗

                        //粒子距離がhより小さい場合だけ計算する
                        if (r2 < h2)
                        {
                            double r = Math.Sqrt(r2); //粒子距離
                            double c = h - r;
                            double n = c / nearParticle.density;
                            sum = Vector3Utility.AddVector3(sum, Vector3Utility.MultiplyScalarVector3(Vector3Utility.SubVector3(nearParticle.velocity, nowParticle.velocity), n));
                        }
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
        Task CalcColiderTerm(List<Particle> particles)
        {
            List<Task> tasks = new();
            for (int i = 0; i < particles.Count; i++)
            {
                int index = i;
                Task task = Task.Run(() =>
                {

                    //reportProgress(reportHeader + "衝突項を計算中..." + i + "/" + particles.length);

                    Particle nowParticle = particles[index]; //今回計算する粒子
                    Vector3 term = new();

                    /*Dictionary<int, bool> filterArray = [];
                    filterPoint(filterArray, data.normalVectors.sort.x, nowParticle.position.x - (-attenuationCoefficient), nowParticle.position.x + (-attenuationCoefficient), "x");
                    filterPoint(filterArray, data.normalVectors.sort.y, nowParticle.position.y - (-attenuationCoefficient), nowParticle.position.y + (-attenuationCoefficient), "y");
                    filterPoint(filterArray, data.normalVectors.sort.z, nowParticle.position.z - (-attenuationCoefficient), nowParticle.position.z + (-attenuationCoefficient), "z");*/
                    for (int j = 0; j < data.Length; j++)
                    {
                        //if (!filterArray[j]) continue;
                        FaceData nowPoint = data[j];
                        if (Is_inside(nowParticle.position, nowPoint))
                        {
                            double distance = Vector3Utility.DotVector3(Vector3Utility.SubVector3(nowParticle.position, nowPoint.CenterOfGravity), nowPoint.NormalVector);
                            Vector3 nowTerm = Vector3Utility.MultiplyScalarVector3(nowPoint.NormalVector, springConstant * distance + attenuationCoefficient * Vector3Utility.DotVector3(nowParticle.velocity, nowPoint.NormalVector));
                            term = Vector3Utility.AddVector3(term, nowTerm);
                        }
                    }
                    terms[index].coliderTerm = term;

                    /*let distance = nowParticle.position.y;
                    if (distance < -attenuationCoefficient) {
                        terms[i].coliderTerm = multiplyScalarVector3(createVector3(0, 1, 0), springConstant * distance + attenuationCoefficient * dotVector3(nowParticle.velocity, createVector3(0, 1, 0)));
                    } else {
                        terms[i].coliderTerm = createVector3();
                    }*/
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
            CalcDensity(_particles).Wait();
            CalcPressure(_particles).Wait();
            CalcPressureTerm(_particles).Wait();
            CalcViscosityTerm(_particles).Wait();
            CalcColiderTerm(_particles).Wait();

            List<Vector3?> tickResult = new();
            for (int i = 0; i < _particles.Count; i++)
            {
                Particle nowParticle = _particles[i];
                //if (nowParticle.is_wall) continue;
                Vector3 a = Vector3Utility.AddVector3(Vector3Utility.AddVector3(Vector3Utility.AddVector3(terms[i].pressureTerm, terms[i].viscosityTerm), terms[i].coliderTerm), g);
                Vector3 v = Vector3Utility.AddVector3(nowParticle.velocity, Vector3Utility.MultiplyScalarVector3(Vector3Utility.AddVector3(nowParticle.acceleration, a), 0.5 * deltaTime));
                Vector3 deltaPosition = Vector3Utility.AddVector3(Vector3Utility.MultiplyScalarVector3(v, deltaTime), Vector3Utility.MultiplyScalarVector3(a, 0.5 * deltaTime * deltaTime));
                nowParticle.acceleration = a;
                nowParticle.velocity = v;
                nowParticle.position = Vector3Utility.AddVector3(nowParticle.position, deltaPosition);
                tickResult.Add(nowParticle.position);
            }
            return tickResult.ToArray();
        }

        void AddParticles(List<Particle> particles)
        {
            Random random = new();
            for (int i = 2; i < 2000; i += 2)
            {
                Particle particle = new();
                particle.position = new Vector3(random.NextDouble() + 0.5, i + 37, 1);
                particle.velocity = new Vector3(0, -10, 0);
                particles.Add(particle);
            }
        }
        public void Start(double simulateSeconds)
        {
            // 初めに実行する処理
            //makeWall(_particles);
            AddParticles(_particles);

            //self.postMessage({ type: "result", content: _particles, time: 0});

            Console.WriteLine("シミュレーション開始");
            for (int i = 0; i < simulateSeconds / deltaTime; i++)
            {
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
                //self.postMessage({ type: "result", content: _particles, time: time});
            }

            //reportProgress("シミュレーション終了");
        }
    }
}
