﻿using System;
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

        internal List<Result> result = [];

        readonly FaceData[] data = data;

        readonly Dictionary<int, Term> terms = [];

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
        readonly static double deltaTime = 0.1/*0.03*/;

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
                        //Console.WriteLine(pressure);
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
            List<Task> tasks = [];
            for (int i = 0; i < particles.Count; i++)
            {
                int index = i;
                Task task = Task.Run(() =>
                {
                    Particle nowParticle = particles[index]; //今回計算する粒子
                    Vector3 term = new();

                    /*Dictionary<int, bool> filterArray = [];
                    filterPoint(filterArray, data.normalVectors.sort.x, nowParticle.position.x - (-dampingCoefficient), nowParticle.position.x + (-dampingCoefficient), "x");
                    filterPoint(filterArray, data.normalVectors.sort.y, nowParticle.position.y - (-dampingCoefficient), nowParticle.position.y + (-dampingCoefficient), "y");
                    filterPoint(filterArray, data.normalVectors.sort.z, nowParticle.position.z - (-dampingCoefficient), nowParticle.position.z + (-dampingCoefficient), "z");*/
                    for (int j = 0; j < data.Length; j++)
                    {
                        //if (!filterArray[j]) continue;
                        FaceData nowPoint = data[j];
                        if (Is_inside(nowParticle.position, nowPoint))
                        {
                            double distance = Math.Abs(Vector3Utility.DotVector3(Vector3Utility.SubVector3(nowParticle.position, nowPoint.CenterOfGravity), nowPoint.NormalVector));
                            Vector3 nowTerm = Vector3Utility.MultiplyScalarVector3(nowPoint.NormalVector, springConstant * distance + dampingCoefficient * Vector3Utility.DotVector3(nowParticle.velocity, nowPoint.NormalVector));
                            term = Vector3Utility.AddVector3(term, nowTerm);
                        }
                    }
                    terms[index].coliderTerm = term;
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

            List<Vector3?> tickResult = [];
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
                //Console.WriteLine(nowParticle.density);
                /*if (terms[i].coliderTerm.Y > 0)
                {
                    Console.WriteLine($"{terms[i].coliderTerm.Y}, {terms[i].pressureTerm.Y}, {terms[i].viscosityTerm.Y}, {_particles[i].acceleration.Y}, {_particles[i].position.Y}");
                }*/
                tickResult.Add(nowParticle.position);
            }
            _particles.RemoveAll(particle => particle.position.Y <= 0);
            return [.. tickResult];
        }

        static void AddParticles(List<Particle> particles)
        {
            for (double z = /*1000*/420; z <= /*1007*/2303; z += particleDistance)
            {
                for (double y = 45; y <= 50; y += particleDistance)
                {
                    Particle particle = new()
                    {
                        position = new Vector3(-700, y, z),
                        velocity = new Vector3(10, 0, 0)
                    };
                    particles.Add(particle);
                }
            }
        }
        public void Start(double simulateSeconds)
        {
            // 初めに実行する処理
            //makeWall(_particles);

            Console.WriteLine("シミュレーション開始");
            for (int i = 0; i < simulateSeconds / deltaTime; i++)
            {
                if (i % 3  == 0)
                {
                    AddParticles(_particles);
                }
                
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
