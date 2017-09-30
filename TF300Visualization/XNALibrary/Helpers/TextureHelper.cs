using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using Microsoft.Xna.Framework.Graphics;

namespace TF300.App.GUI.DatabaseUI.XNALibrary.Helpers
{
    using Color = Microsoft.Xna.Framework.Graphics.Color;
    public class TextureHelper
    {

        public static Texture2D Texture2DFromBitmap(Bitmap bmp, GraphicsDevice graphicsDevice)
        {
            Color[] pixels = new Color[bmp.Width * bmp.Height];
            for (int y = 0; y < bmp.Height; y++)
            {
                for (int x = 0; x < bmp.Width; x++)
                {
                    System.Drawing.Color c = bmp.GetPixel(x, y);
                    pixels[(y * bmp.Width) + x] = new Color(c.R, c.G, c.B, c.A);
                }
            }

            Texture2D myTex = new Texture2D(
              graphicsDevice,
              bmp.Width,
              bmp.Height,
              1,
              TextureUsage.None,
              SurfaceFormat.Color);

            myTex.SetData<Color>(pixels);
            return myTex;
        }
    }

}
