using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using XNAHelper;
using XNAHelper.Interpolaters;
using System.Collections.Generic;
using XNAHelper.CustomVertexs;
using XNAHelper.Helpers;
using TF300.App.GUI.DatabaseUI;
using System;
using Microsoft.Xna.Framework.Content;
using System.Drawing;
using TF300.App.GUI.DatabaseUI.XNALibrary.Helpers;


namespace DrawUserPrimitives.Primitives
{
    using ColorD = System.Drawing.Color;
    using Color = Microsoft.Xna.Framework.Graphics.Color;
    using FillMode = Microsoft.Xna.Framework.Graphics.FillMode;
    using System.Drawing.Drawing2D;

    /// <summary>
    /// Geometric primitive class for drawing cubes.
    /// </summary>
    public class DrawerColorMap : Drawer, IDisposable
    {
        bool _Flat;

        public bool Flat
        {
            get { return _Flat; }
            set
            {
                _Flat = value;
                _OptionApplyed = false;
            }
        }

        VertexPositionColorNormal[] mVertices;
        ushort[] mIndices;

        VertexPositionColorNormal[] mDrawVertices;


        VertexBuffer vertexBufferBody;
        IndexBuffer indexBufferBody;


        float mMinValue;
        float mMaxValue;
        //Matrix matrixScale;
        /// <summary>
        /// Constructs a new cube primitive, with the specified size.
        /// </summary>
        public DrawerColorMap(GraphicsDevice graphicsDevice, ContentManager Content, VertexPositionColorNormal[] verticies, ushort[] indices)
            : base(graphicsDevice, Content)
        {

            //图例
            InitializeLegend();

            //顶点处理
            mVertices = new VertexPositionColorNormal[verticies.Length];
            verticies.CopyTo(mVertices, 0);
            mIndices = new ushort[indices.Length];
            indices.CopyTo(mIndices, 0);

            PrimitiveHelper.RenderGradientColor(mVertices, Color.Yellow, Color.Red, Color.Green);



            float minValue = float.MaxValue;
            float maxValue = float.MinValue;
            foreach (VertexPositionColorNormal Vertex in verticies)
            {
                if (minValue > Vertex.Position.Z) minValue = Vertex.Position.Z;
                if (maxValue < Vertex.Position.Z) maxValue = Vertex.Position.Z;
            }
            mMinValue = minValue;
            mMaxValue = maxValue;

            PrimitiveHelper.ScaleZToRange(mVertices, 10f, minValue, maxValue);

            //在缩放了Z之后计算法线才可得到正确的光照

            //int[] indices = MeshGenerator.CreateCircleIndex(50, 50);
            PrimitiveHelper.GenerateNormalsForTriangleList(mVertices, indices);


            //InitializeEffect();
        }

        Texture2D mLegendTexture;
        int mLegendWidth = 30;
        int mLegendHeight = 400;

        private void InitializeLegend()
        {
            //生成图例的贴图纹理

            System.Drawing.Bitmap bmp = new System.Drawing.Bitmap(
                                               mLegendWidth, mLegendHeight,
                                               System.Drawing.Imaging.PixelFormat.Format24bppRgb
                                             );
            Graphics graphics = Graphics.FromImage(bmp);
            LinearGradientBrush brush = new LinearGradientBrush(new System.Drawing.Rectangle(0, 0, mLegendWidth, mLegendHeight/2), ColorD.Red, ColorD.Yellow, LinearGradientMode.Vertical);
            graphics.FillRectangle(brush, new System.Drawing.Rectangle(0, 0, mLegendWidth, mLegendHeight / 2));
            ColorD[] colors = new ColorD[2];
            colors[0] = ColorD.Yellow;
            colors[1] = ColorD.Green;
            brush.LinearColors = colors;
            graphics.FillRectangle(brush, new System.Drawing.Rectangle(0, mLegendHeight / 2, mLegendWidth, mLegendHeight / 2));
            brush.Dispose();
            graphics.Dispose();
            //位图转换到纹理

            mLegendTexture = TextureHelper.Texture2DFromBitmap(bmp, this.mGraphicsDevice);

        }

        private Vector2 mLegendPosition = new Vector2(20,120);

        /// <summary>
        /// 绘制图型
        /// </summary>
        protected override void Draw(BasicEffect effect)
        {

            GraphicsDevice graphicsDevice = effect.GraphicsDevice;
            RenderState renderState = graphicsDevice.RenderState;
            renderState.FillMode = FillMode.Solid;
            renderState.CullMode = CullMode.None;

            //设置顶点与索引缓存

            graphicsDevice.VertexDeclaration = vertexDeclaration;
            //graphicsDevice.Vertices[0].SetSource(vertexBufferContour, 0, VertexPositionColorNormal.SizeInBytes);
            //graphicsDevice.Indices = indexBufferContour;


            effect.Begin();

            foreach (EffectPass effectPass in effect.CurrentTechnique.Passes)
            {
                effectPass.Begin();

                graphicsDevice.Vertices[0].SetSource(vertexBufferBody, 0, VertexPositionColorNormal.SizeInBytes);
                graphicsDevice.Indices = indexBufferBody;
                graphicsDevice.DrawIndexedPrimitives(PrimitiveType.TriangleList, 0, 0, this.mDrawVertices.Length, 0, mIndices.Length / 3);
                effectPass.End();
            }

            mSpriteBatch.Begin();
            mSpriteBatch.Draw(mLegendTexture, mLegendPosition, Color.White);
            mSpriteBatch.DrawString(mFont, mMaxValue.ToString("f5"), new Vector2(mLegendPosition.X + mLegendWidth, mLegendPosition.Y), Color.SteelBlue, 0, new Vector2(0, 0), 0.7f, SpriteEffects.None, 0);
            mSpriteBatch.DrawString(mFont, mMinValue.ToString("f5"), new Vector2(mLegendPosition.X + mLegendWidth, mLegendPosition.Y + mLegendHeight), Color.SteelBlue, 0, new Vector2(0, 0), 0.7f, SpriteEffects.None, 0);
            mSpriteBatch.DrawString(mFont, ((mMaxValue+mMinValue)/2).ToString("f5"), new Vector2(mLegendPosition.X + mLegendWidth, mLegendPosition.Y + mLegendHeight/2), Color.SteelBlue, 0, new Vector2(0, 0), 0.7f, SpriteEffects.None, 0);
            mSpriteBatch.End();


            effect.End();
        }

        protected override void ApplyOptions()
        {
            ClearBuffer();
            mDrawVertices = new VertexPositionColorNormal[mVertices.Length];
            mVertices.CopyTo(mDrawVertices, 0);

            if (_Flat)
            {
                for (int i = 0; i < mDrawVertices.Length; i++)
                {
                    mDrawVertices[i].Position.Z = 0;
                    mDrawVertices[i].Normal = Vector3.Up;
                }
            }
            else
            {
            }

            vertexBufferBody = new VertexBuffer(mGraphicsDevice, typeof(VertexPositionColorNormal), this.mDrawVertices.Length, BufferUsage.None);
            vertexBufferBody.SetData(mDrawVertices);
            indexBufferBody = new IndexBuffer(mGraphicsDevice, typeof(ushort), this.mIndices.Length, BufferUsage.None);
            indexBufferBody.SetData(mIndices);

            vertexDeclaration = new VertexDeclaration(mGraphicsDevice, VertexPositionColorNormal.VertexElements);
            _OptionApplyed = true;


        }

        protected override void ClearBuffer()
        {
            if (vertexBufferBody != null)
            {
                vertexBufferBody.Dispose();
            }
            if (indexBufferBody != null)
            {
                indexBufferBody.Dispose();
            }

            base.ClearBuffer();
        }

    }
}
