using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using DisasterSimulation;
using static System.Runtime.InteropServices.JavaScript.JSType;
using static DisasterSimulation.Simulator;

namespace DisasterSimulation
{
    internal class FaceData
    {
        public int index;
        public Vector3[] triangle = new Vector3[3];
        public Vector3 centerOfGravity = new();
        public Vector3 normalVector = new();
        public List<InsideJudge> insideJudge = new();
        public class InsideJudge
        {
            public Vector3 point = new();
            public Vector3 normalVector = new();
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

        readonly FaceData[] data = data;

        Dictionary<int, Term> terms = new();

        readonly static double h = 0.3f/* 0.012f */; //影響半径
        readonly static double particleMass = 0.0002f; //粒子の質量

        readonly static Vector3 g = new(0, -9.8f, 0);  // 重力加速度
        readonly static double pressureStiffness = 200; //圧力係数
        readonly static double restDensity = 1000; //静止密度
        readonly static double viscosity = 0.000001f;  // 粘性係数
        readonly static double attenuationCoefficient = -5;  // ダンパ係数（javascript側と同じ値を用いること）
        readonly static double springConstant = -5;  // ばね係数

        readonly double densityCoef = particleMass * 315 / (64 * Math.PI * Math.Pow(h, 9)); //密度計算で使う

        readonly double pressureCoef = particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //圧力項計算で使う
        readonly double viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.Pow(h, 6)); //粘性項計算で使う

        List<Particle> _particles = new();
        readonly static double deltaTime = 0.01;

        void CalcDensity(List<Particle> particles)
        {
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算

                //reportProgress(reportHeader + "密度を計算中..." + i + "/" + particles.Count);

                Particle nowParticle = particles[i]; //今回計算する粒子
                double sum = 0; //足し合わせる変数
                for (int j = 0; j < particles.Count; j++)
                { //他の粒子全てについて
                    if (i == j) { continue; } //自分自身だったらスキップ
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
            }
        }

        void CalcPressure(List<Particle> particles)
        {
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の圧力を計算

                //reportProgress(reportHeader + "圧力を計算中..." + i + "/" + particles.Count);

                particles[i].pressure = pressureStiffness * (particles[i].density - restDensity);
            }
        }

        void CalcPressureTerm(List<Particle> particles)
        {
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算

                //reportProgress(reportHeader + "圧力項を計算中..." + i + "/" + particles.length);

                Particle nowParticle = particles[i]; //今回計算する粒子
                //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                Vector3 sum = new(); //足し合わせる変数
                for (int j = 0; j < particles.Count; j++)
                { //他の粒子全てについて
                    if (i == j) { continue; } //自分自身だったらスキップ
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

                if (!terms.ContainsKey(i)) terms[i] = new();
                terms[i].pressureTerm = Vector3Utility.MultiplyScalarVector3(sum, (-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
            }
        }

        void CalcViscosityTerm(List<Particle> particles)
        {
            double h2 = h * h; //事前にhの二乗を計算しておく
            for (int i = 0; i < particles.Count; i++)
            { //一つづつ粒子の密度を計算

                //reportProgress(reportHeader + "粘性項を計算中..." + i + "/" + particles.length);

                Particle nowParticle = particles[i]; //今回計算する粒子
                //if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
                Vector3 sum = new(); //足し合わせる変数
                for (int j = 0; j < particles.Count; j++)
                { //他の粒子全てについて
                    if (i == j) { continue; } //自分自身だったらスキップ
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

                if (terms.ContainsKey(i)) terms[i] = new();
                terms[i].viscosityTerm = Vector3Utility.MultiplyScalarVector3(sum, viscosityCoef);  // 粘性項が求まった
            }
        }

        static bool Is_inside(Vector3 position, FaceData point)
        {
            for (int i = 0; i < point.insideJudge.Count; i++)
            {
                double dot = Vector3Utility.DotVector3(point.insideJudge[i].normalVector, Vector3Utility.SubVector3(position, point.insideJudge[i].point));
                if (dot < 0) return false;
            }
            return true;
            /*const dot = dotVector3(point.normalVector, subVector3(position, point.centerOfGravity));
            return dot <= 0;*/
        }
        void CalcColiderTerm(List<Particle> particles)
        {
            for (int i = 0; i < particles.Count; i++)
            {

                //reportProgress(reportHeader + "衝突項を計算中..." + i + "/" + particles.length);

                Particle nowParticle = particles[i]; //今回計算する粒子
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
                        double distance = Vector3Utility.DotVector3(Vector3Utility.SubVector3(nowParticle.position, nowPoint.centerOfGravity), nowPoint.normalVector);
                        Vector3 nowTerm = Vector3Utility.MultiplyScalarVector3(nowPoint.normalVector, springConstant * distance + attenuationCoefficient * Vector3Utility.DotVector3(nowParticle.velocity, nowPoint.normalVector));
                        term = Vector3Utility.AddVector3(term, nowTerm);
                    }
                }
                terms[i].coliderTerm = term;

                /*let distance = nowParticle.position.y;
                if (distance < -attenuationCoefficient) {
                    terms[i].coliderTerm = multiplyScalarVector3(createVector3(0, 1, 0), springConstant * distance + attenuationCoefficient * dotVector3(nowParticle.velocity, createVector3(0, 1, 0)));
                } else {
                    terms[i].coliderTerm = createVector3();
                }*/
            }
        }

        void Tick()
        {
            CalcDensity(_particles);
            CalcPressure(_particles);
            CalcPressureTerm(_particles);
            CalcViscosityTerm(_particles);
            CalcColiderTerm(_particles);

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
            }
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

            for (int i = 0; i < simulateSeconds / deltaTime; i++)
            {
                //double time = (i + 1) * deltaTime;
                //reportHeader = time + "秒目/" + simulateSeconds + "秒 ";
                Tick();
                //self.postMessage({ type: "result", content: _particles, time: time});
            }

            //reportProgress("シミュレーション終了");
        }
    }
}
