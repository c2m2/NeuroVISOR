using System;
using Vector = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MathNet.Numerics.LinearAlgebra;

namespace C2M2.NeuronalDynamics.Simulation
{
    public static class PremadeChannels
    {
        /// <summary>
        /// Potassium channel matching the solver's definitions verbatim
        /// gK = 5.0e1, eK = -90e-3
        /// alpha_n, beta_n exactly as in SparseSolverTestv1
        /// </summary>
        public static IonChannel PotassiumChannel(int nodeCount)
        {
            double gk = 5.0 * 1.0E1;   
            double ek = -90.0 * 1.0E-3;

            IonChannel potassiumChannel = new IonChannel("Potassium Channel", gk, ek);

            Func<Vector, Vector> alpha_n = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin); 
                return (1.0E3) * (0.032)
                    * (15.0 - Vin)
                    .PointwiseDivide(((15.0 - Vin) / 5.0).PointwiseExp() - 1.0);
            };

            Func<Vector, Vector> beta_n = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.5)
                    * ((10.0 - Vin) / 40.0)
                    .PointwiseExp();
            };

            potassiumChannel.AddGatingVariable(
                new GatingVariable("n", alpha_n, beta_n, 4, 0.0376969, nodeCount)
            );

            return potassiumChannel;
        }

        /// <summary>
        /// Sodium channel matching the solver's definitions verbatim
        /// gNa = 60.0e1, eNa = 50.0e-3
        /// alpha_m, beta_m, alpha_h, beta_h exactly as in SparseSolverTestv1
        /// </summary>
        public static IonChannel SodiumChannel(int nodeCount)
        {
            double gna = 60.0 * 1.0E1;    // 600 S/m²
            double ena = 50.0 * 1.0E-3;   // 0.050 V

            IonChannel sodiumChannel = new IonChannel("Sodium Channel", gna, ena);

            // alpha_m(V)
            Func<Vector, Vector> alpha_m = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.32)
                    * (13.0 - Vin)
                    .PointwiseDivide(((13.0 - Vin) / 4.0).PointwiseExp() - 1.0);
            };

            Func<Vector, Vector> beta_m = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.28)
                    * (Vin - 40.0)
                    .PointwiseDivide(((Vin - 40.0) / 5.0).PointwiseExp() - 1.0);
            };

            Func<Vector, Vector> alpha_h = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.128)
                    * ((17.0 - Vin) / 18.0)
                    .PointwiseExp();
            };

            Func<Vector, Vector> beta_h = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * 4.0
                    / (((40.0 - Vin) / 5.0).PointwiseExp() + 1.0);
            };

            sodiumChannel.AddGatingVariable(
                new GatingVariable("m", alpha_m, beta_m, 3, 0.0147567, nodeCount)
            );
            sodiumChannel.AddGatingVariable(
                new GatingVariable("h", alpha_h, beta_h, 1, 0.9959410, nodeCount)
            );

            return sodiumChannel;
        }

        /// <summary>
        /// Calcium channel matching the solver's definitions verbatim
        /// gCa = 1.0e1, eCa = 120.0e-3
        /// alpha_q, beta_q, alpha_r, beta_r exactly as in SparseSolverTestv1
        /// </summary>
        public static IonChannel CalciumChannel(int nodeCount)
        {
            double gca = 1.0 * 1.0E1;
            double eca = 120.0 * 1.0E-3;

            IonChannel calciumChannel = new IonChannel("Calcium Channel", gca, eca);

            Func<Vector, Vector> alpha_q = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.055)
                    * (27.0 - Vin)
                    / (((-27.0 - Vin) / 3.8).PointwiseExp() - 1.0);
            };

            Func<Vector, Vector> beta_q = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.94)
                    * (((-75.0 - Vin) / 17.0).PointwiseExp());
            };

            Func<Vector, Vector> alpha_r = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.000457)
                    * (((-13.0 - Vin) / 50.0).PointwiseExp());
            };

            Func<Vector, Vector> beta_r = voltage =>
            {
                var Vin = voltage.Clone();
                Vin.Multiply(1.0E3, Vin);
                return (1.0E3) * (0.0065)
                    / (((-15.0 - Vin) / 28.0).PointwiseExp() + 1.0);
            };

            calciumChannel.AddGatingVariable(
                new GatingVariable("q", alpha_q, beta_q, 2, 0.0, nodeCount)
            );
            calciumChannel.AddGatingVariable(
                new GatingVariable("r", alpha_r, beta_r, 1, 0.0, nodeCount)
            );

            return calciumChannel;
        }

        /// <summary>
        /// Chloride channel matching the solver's definitions verbatim
        /// user can pass gCl, eCl
        /// alpha_cl, beta_cl as in SparseSolverTestv1
        /// </summary>
        // public static IonChannel ChlorideChannel(int nodeCount, double gcl, double ecl)
        // {
        //     IonChannel chlorideChannel = new IonChannel("Chloride Channel", gcl, ecl);

        //     Func<Vector, Vector> alpha_cl = voltage =>
        //     {
        //         var Vin = voltage.Clone();
        //         Vin.Multiply(1.0E3, Vin);
        //         return (1.0E3 * 0.07)
        //             * (Vin.Add(20.0))
        //             .PointwiseDivide(
        //                 (Vin.Add(20.0).Divide(10.0))
        //                 .PointwiseExp()
        //                 .Subtract(1.0)
        //             );
        //     };

        //     Func<Vector, Vector> beta_cl = voltage =>
        //     {
        //         var Vin = voltage.Clone();
        //         Vin.Multiply(1.0E3, Vin);
        //         return (1.0E3 * 0.1)
        //             * (-(Vin.Add(30.0)).Divide(10.0))
        //             .PointwiseExp();
        //     };

        //     chlorideChannel.AddGatingVariable(
        //         new GatingVariable("x", alpha_cl, beta_cl, 1, 0.5, nodeCount)
        //     );

        //     return chlorideChannel;
        // }

        // // /// <summary>
        // // /// Leakage channel matching the solver's gl=0.0, el=-70e-3
        // // /// </summary>
        
        public static IonChannel LeakageChannel(int nodeCount)
        {
            double gl = 0.0;
            double el = -70.0 * 1.0E-3;
            return new IonChannel("Leakage Channel", gl, el);
        }
        
    }
}
