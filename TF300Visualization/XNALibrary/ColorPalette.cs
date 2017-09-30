using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework;

namespace TF300.App.GUI.DatabaseUI.XNALibrary
{
    public class ColorPalette
    {
        List<Vector3> mListColorKey = new List<Vector3>();

        public ColorPalette()
        {

            mListColorKey.Add(new Vector3(0, 0, 127));
            mListColorKey.Add(new Vector3(0, 0, 255));
            mListColorKey.Add(new Vector3(0, 165, 255));
            mListColorKey.Add(new Vector3(0, 255, 255));
            mListColorKey.Add(new Vector3(0, 255, 165));
            mListColorKey.Add(new Vector3(0, 255, 0));
            mListColorKey.Add(new Vector3(165, 255, 0));
            mListColorKey.Add(new Vector3(255, 255, 0));
            mListColorKey.Add(new Vector3(255, 165, 0));
            mListColorKey.Add(new Vector3(255, 0, 0));
            mListColorKey.Add(new Vector3(127, 0, 0));


        }

        public Color GetColor(float amount)
        {
            int nStep = (int)(amount / (1f / (mListColorKey.Count-1)));
            if (amount == 1f) nStep -= 1;
            float stepAmount = amount - (1f / (mListColorKey.Count-1)) * nStep;

            Vector3 v0 = mListColorKey[nStep];
            Vector3 v1 = mListColorKey[nStep+1];

            Vector3 vColor = Vector3.Lerp(v0, v1, stepAmount);

            Color color = new Color((byte)vColor.X,(byte)vColor.Y,(byte)vColor.Z);
            
            return color;
            //return new Color( (vColor.X,vColor.Y,vColor.Z);

        }
    }
}
