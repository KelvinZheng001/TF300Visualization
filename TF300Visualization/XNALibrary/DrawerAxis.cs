using System;
using System.Collections.Generic;
using System.Text;
using XNAHelper.CustomVertexs;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using XNAHelper.Interpolaters;
using Microsoft.Xna.Framework.Content;
using XNAHelper.Helpers;

namespace TF300.App.GUI.DatabaseUI.XNALibrary
{
    public class DrawerAxis : Drawer
    {

        bool _Flat;
        /// <summary>
        /// 指定是否绘制为平面
        /// </summary>
        public bool Flat
        {
            get { return _Flat; }
            set
            {
                _Flat = value;
                _OptionApplyed = false;
            }
        }
        private VertexPositionColorNormal[] mVerticies = new VertexPositionColorNormal[4];
        private short[] mIndices = new short[8];
        private VertexPositionColorNormal[] mIndicatorVerticies = new VertexPositionColorNormal[3];
        

        private VertexPositionColorNormal[] mDataPointVerticies;
        private VertexPositionColorNormal[] mDrawDataPointVerticies;


        public DrawerAxis(GraphicsDevice graphicsDevice, ContentManager Content, VertexPositionColorNormal[] dataVerticies, VertexPositionColorNormal[] verticies)
            : base(graphicsDevice, Content)
        {
            mVerticies[0] = new VertexPositionColorNormal(new Vector3(-25f, -25f, 0f), Color.Black, Vector3.Up);
            mVerticies[1] = new VertexPositionColorNormal(new Vector3(-25f, 25f, 0f), Color.Black, Vector3.Up);
            mVerticies[2] = new VertexPositionColorNormal(new Vector3(25f, 25f, 0f), Color.Black, Vector3.Up);
            mVerticies[3] = new VertexPositionColorNormal(new Vector3(25f, -25f, 0f), Color.Black, Vector3.Up);
            mIndices[0] = 0;
            mIndices[1] = 1;
            mIndices[2] = 1;
            mIndices[3] = 2;
            mIndices[4] = 2;
            mIndices[5] = 3;
            mIndices[6] = 3;
            mIndices[7] = 0;

            mIndicatorVerticies[0] = new VertexPositionColorNormal(new Vector3(0f, -25f, 0f), Color.Black, Vector3.Up);
            mIndicatorVerticies[1] = new VertexPositionColorNormal(new Vector3(-1f, -27f, 0f), Color.Black, Vector3.Up);
            mIndicatorVerticies[2] = new VertexPositionColorNormal(new Vector3(1f, -27f, 0f), Color.Black, Vector3.Up);

            
            mDataPointVerticies = new VertexPositionColorNormal[dataVerticies.Length];
            dataVerticies.CopyTo(mDataPointVerticies, 0);

            //计算最大最小高度
            float minValue = float.MaxValue;
            float maxValue = float.MinValue;
            foreach (VertexPositionColorNormal Vertex in verticies)
            {
                if (minValue > Vertex.Position.Z) minValue = Vertex.Position.Z;
                if (maxValue < Vertex.Position.Z) maxValue = Vertex.Position.Z;
            }


            PrimitiveHelper.ScaleZToRange(mDataPointVerticies, 10f, minValue, maxValue);

        }

        protected override void ApplyOptions()
        {
            vertexDeclaration = new VertexDeclaration(mGraphicsDevice, VertexPositionColorNormal.VertexElements);
            mDrawDataPointVerticies = new VertexPositionColorNormal[mDataPointVerticies.Length];
            mDataPointVerticies.CopyTo(mDrawDataPointVerticies, 0);
            if (_Flat == true)
            {
                for (int i = 0; i < mDrawDataPointVerticies.Length; i++)
                {
                    mDrawDataPointVerticies[i].Position.Z = 0f;
                }
            }
        }
        float _Yaw;
        float _Roll;
        float _Pitch;

        /// <summary>
        /// 设置旋转角度的信息
        /// </summary>
        public void SetYawRollPitchInfo(float yaw, float roll, float pitch)
        {
            _Yaw = yaw;
            _Roll = roll;
            _Pitch = pitch;
        }


        protected override void Draw(BasicEffect effect)
        {
            GraphicsDevice graphicsDevice = effect.GraphicsDevice;
            RenderState renderState = graphicsDevice.RenderState;
            renderState.FillMode = FillMode.Solid;
            renderState.CullMode = CullMode.None;



            //设置顶点与索引缓存
            graphicsDevice.VertexDeclaration = vertexDeclaration;
            effect.Begin();

            //绘制数据点十字标记
            VertexPositionColorNormal[] verticiesLink = new VertexPositionColorNormal[mDrawDataPointVerticies.Length * 4];
            Vector3 location;
            //Vector3 point1 = graphicsDevice.Viewport.Unproject(new Vector3(500f, 500f, 0f), effect.Projection, effect.View, effect.World);
            for (int i = 0; i < mDrawDataPointVerticies.Length; i++)
            {
                location = graphicsDevice.Viewport.Project(mDrawDataPointVerticies[i].Position, effect.Projection, effect.View, effect.World);
                Vector3 point1 = graphicsDevice.Viewport.Unproject(new Vector3(location.X - 5, location.Y, 0f), effect.Projection, effect.View, effect.World);
                Vector3 point2 = graphicsDevice.Viewport.Unproject(new Vector3(location.X + 5, location.Y, 0f), effect.Projection, effect.View, effect.World);
                Vector3 point3 = graphicsDevice.Viewport.Unproject(new Vector3(location.X, location.Y - 5, 0f), effect.Projection, effect.View, effect.World);
                Vector3 point4 = graphicsDevice.Viewport.Unproject(new Vector3(location.X, location.Y + 5, 0f), effect.Projection, effect.View, effect.World);

                verticiesLink[i * 4] = new VertexPositionColorNormal(point1, Color.YellowGreen, Vector3.Up);
                verticiesLink[i * 4 + 1] = new VertexPositionColorNormal(point2, Color.YellowGreen, Vector3.Up);
                verticiesLink[i * 4 + 2] = new VertexPositionColorNormal(point3, Color.YellowGreen, Vector3.Up);
                verticiesLink[i * 4 + 3] = new VertexPositionColorNormal(point4, Color.YellowGreen, Vector3.Up);
            }
            short[] indexLink = new short[verticiesLink.Length];
            for (short i = 0; i < indexLink.Length; i++)
            {
                indexLink[i] = i;
            }


            foreach (EffectPass effectPass in effect.CurrentTechnique.Passes)
            {
                effectPass.Begin();
                if (_Flat == false)
                    mGraphicsDevice.DrawUserIndexedPrimitives<VertexPositionColorNormal>(PrimitiveType.LineList, mVerticies, 0, mVerticies.Length, mIndices, 0, mIndices.Length / 2);
                mGraphicsDevice.DrawUserPrimitives<VertexPositionColorNormal>(PrimitiveType.TriangleList, mIndicatorVerticies, 0, 1);
                //连线
                mGraphicsDevice.DrawUserIndexedPrimitives<VertexPositionColorNormal>(PrimitiveType.LineList, verticiesLink, 0, verticiesLink.Length, indexLink, 0, indexLink.Length / 2);

                effectPass.End();
            }
            effect.End();



            //绘制指示箭头与测量点文字

            _Yaw = MathHelper.ToDegrees(_Yaw);
            _Roll = MathHelper.ToDegrees(_Roll);
            _Pitch = MathHelper.ToDegrees(_Pitch);
            string strYaw = string.Format("Yaw:{0:f3}",_Yaw);
            string strRoll = string.Format("Roll:{0:f3}", _Roll);
            string strPitch = string.Format("Pitch:{0:f3}", _Pitch);


            mSpriteBatch.Begin();
            location = graphicsDevice.Viewport.Project(mIndicatorVerticies[2].Position, effect.Projection, effect.View, effect.World);
            mSpriteBatch.DrawString(mFont, "Indicator", new Vector2(location.X, location.Y), Color.Blue);
            for (int i = 0; i < mDrawDataPointVerticies.Length; i++)
            {
                location = graphicsDevice.Viewport.Project(mDrawDataPointVerticies[i].Position, effect.Projection, effect.View, effect.World);
                mSpriteBatch.DrawString(mFont, string.Format("P{0}", i, mDrawDataPointVerticies[i].Position.Z), new Vector2(location.X + 7, location.Y), Color.DarkGray);
            }
            try
            {
                mSpriteBatch.DrawString(mFont, strRoll, new Vector2(20, 40), Color.SteelBlue);
                mSpriteBatch.DrawString(mFont, strYaw, new Vector2(20, 60), Color.SteelBlue);
                mSpriteBatch.DrawString(mFont, strPitch, new Vector2(20, 80), Color.SteelBlue);
            }
            catch
            {
            }
            mSpriteBatch.End();











        }
    }
}
