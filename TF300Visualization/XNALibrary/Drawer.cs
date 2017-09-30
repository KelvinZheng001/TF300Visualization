using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Content;


namespace TF300.App.GUI.DatabaseUI
{
    public class Drawer: IDisposable
    {

        
        protected BasicEffect basicEffect;
        protected VertexDeclaration vertexDeclaration;
        protected GraphicsDevice mGraphicsDevice = null;
        protected bool _OptionApplyed = false;
        protected ContentManager mContentManager;
        protected SpriteBatch mSpriteBatch;
        protected SpriteFont mFont;

        public Drawer(GraphicsDevice graphicsDevice, ContentManager Content)
        {
            mContentManager = Content;
            InitializeEffect(graphicsDevice);
            
        }

        /// <summary>
        /// 初始化绘图与顶点缓冲
        /// </summary>
        /// <param name="graphicsDevice"></param>
        protected void InitializeEffect(GraphicsDevice graphicsDevice)
        {

            // Create a BasicEffect, which will be used to render the primitive.
            mGraphicsDevice = graphicsDevice;
            basicEffect = new BasicEffect(graphicsDevice, null);
            //basicEffect.EnableDefaultLighting();
            //basicEffect.LightingEnabled = false;
            basicEffect.PreferPerPixelLighting = true;
            mSpriteBatch = new SpriteBatch(mGraphicsDevice);
            mFont = mContentManager.Load<SpriteFont>("Arial");
        }

        /// <summary>
        /// 设置旋转,效果
        /// </summary>
        public void Draw(Matrix world, Matrix view, Matrix projection)
        {
            if (_OptionApplyed == false)
            {
                ApplyOptions();
            }


            // Set BasicEffect parameters.
            basicEffect.World = world;
            basicEffect.View = view;
            basicEffect.Projection = projection;
            basicEffect.VertexColorEnabled = true;
            //basicEffect.DiffuseColor = color.ToVector3();
            //basicEffect.Alpha = color.A / 255.0f;

            // Set important renderstates.
            RenderState renderState = basicEffect.GraphicsDevice.RenderState;

            renderState.AlphaTestEnable = false;
            renderState.DepthBufferEnable = true;
            renderState.DepthBufferFunction = CompareFunction.LessEqual;

            //if (color.A < 255)
            //{
            //    // Set renderstates for alpha blended rendering.
            //    renderState.AlphaBlendEnable = true;
            //    renderState.AlphaBlendOperation = BlendFunction.Add;
            //    renderState.SourceBlend = Blend.SourceAlpha;
            //    renderState.DestinationBlend = Blend.InverseSourceAlpha;
            //    renderState.SeparateAlphaBlendEnabled = false;
            //    renderState.DepthBufferWriteEnable = false;
            //}
            //else
            //{
            // Set renderstates for opaque rendering.
            renderState.AlphaBlendEnable = false;
            renderState.DepthBufferWriteEnable = true;
            //}

            // Draw the model, using BasicEffect.
            Draw(basicEffect);
        }

        protected virtual void Draw(BasicEffect effect)
        {
        }

        protected virtual void ApplyOptions()
        {
        }

        protected virtual void ClearBuffer()
        {
            if (vertexDeclaration != null)
                vertexDeclaration.Dispose();

        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposing)
            {
                ClearBuffer();
                if (basicEffect != null)
                    basicEffect.Dispose();
            }
        }
    }
}
