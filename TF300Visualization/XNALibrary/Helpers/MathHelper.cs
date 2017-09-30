using System;
using System.Collections.Generic;
using System.Text;
//using TF300.App.GUI.DatabaseUI.Utility;
using Microsoft.Xna.Framework;

namespace TF300.App.GUI.DatabaseUI.XNALibrary.Helpers
{
    public class MathHelper
    {
        public static void ExtractYarRollPithFromMatrix(Microsoft.Xna.Framework.Matrix matrix, ref float yaw, ref float roll, ref float pitch)
        {

            Vector3 X, Y, Z;
            X.X = matrix.M11;
            X.Y = matrix.M21;
            X.Z = matrix.M31;
            Y.X = matrix.M12;
            Y.Y = matrix.M22;
            Y.Z = matrix.M32;
            Z.X = matrix.M13;
            Z.Y = matrix.M23;
            Z.Z = matrix.M33;

            
            
            // get angles from basis
            pitch = (float)Math.Atan2(Z.Y , Z.Z);
            yaw = (float)Math.Asin(-Z.X);
            roll = (float)Math.Atan2(Y.X , X.X);

        }

    }
}
