using System;
using System.Collections.Generic;
using System.Text;
using WinFormsGraphicsDevice;
using System.Windows.Forms;
using TF300.App.GUI.DatabaseUI.XNALibrary;
using XNAHelper;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework;
using DrawUserPrimitives.Primitives;
using System.Diagnostics;
using XNAHelper.Interpolaters;
using System.Drawing;
using System.IO;

namespace TF300Visualization
{
    using Color = Microsoft.Xna.Framework.Graphics.Color;
    using DrawUserPrimitives;
    using XNAHelper.CustomVertexs;
    using XNAHelper.Helpers;
    using TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters;

    public class XNAMesh3D : GraphicsDeviceControl
    {

        public event EventHandler ImageCaptured;
        GraphicTypeEnum _GrpahicType;

        public GraphicTypeEnum GrpahicType
        {
            get { return _GrpahicType; }
            set
            {
                _GrpahicType = value;
                PrepareGraphic();

            }
        }


        int _CountourLayerCount = 10;

        public int CountourLayerCount
        {
            get { return _CountourLayerCount; }
            set { _CountourLayerCount = value; }
        }

        public void SetCountourLayerCount(int nLayerCount)
        {
            _CountourLayerCount = nLayerCount;
            mDrawerContour = null;
            PrepareGraphic();
        }


        private void PrepareGraphic()
        {
            if (mDrawerAxis == null)
                mDrawerAxis = new DrawerAxis(GraphicsDevice, Content, mDataPointVertices, mVertices);
            switch (_GrpahicType)
            {
                case GraphicTypeEnum.ContourLine2D:
                    if (mDrawerContour == null)
                    {
                        mDrawerContour = new DrawerContourMap(this.GraphicsDevice, Content, mVertices, mIndices, mEdgeVertices, CountourLayerCount);
                    }
                    mDrawerContour.Flat = true;
                    mDrawerAxis.Flat = true;
                    mDrawerContour.DrawContourLine = true;
                    mDrawerContour.DrawContourZone = false;
                    CenterView();
                    break;
                case GraphicTypeEnum.ContourLine3D:
                    if (mDrawerContour == null)
                    {
                        mDrawerContour = new DrawerContourMap(this.GraphicsDevice, Content, mVertices, mIndices, mEdgeVertices, CountourLayerCount);
                    }
                    mDrawerContour.Flat = false;
                    mDrawerAxis.Flat = false;
                    mDrawerContour.DrawContourLine = true;
                    mDrawerContour.DrawContourZone = false;
                    break;
                case GraphicTypeEnum.ContourZone2D:
                    if (mDrawerContour == null)
                    {
                        mDrawerContour = new DrawerContourMap(this.GraphicsDevice, Content, mVertices, mIndices, mEdgeVertices, CountourLayerCount);
                    }
                    mDrawerContour.Flat = true;
                    mDrawerAxis.Flat = true;
                    mDrawerContour.DrawContourLine = false;
                    mDrawerContour.DrawContourZone = true;
                    CenterView();
                    break;
                case GraphicTypeEnum.ContourZone3D:
                    if (mDrawerContour == null)
                    {
                        mDrawerContour = new DrawerContourMap(this.GraphicsDevice, Content, mVertices, mIndices, mEdgeVertices, CountourLayerCount);
                    }
                    mDrawerContour.Flat = false;
                    mDrawerAxis.Flat = false;
                    mDrawerContour.DrawContourLine = false;
                    mDrawerContour.DrawContourZone = true;
                    break;
                case GraphicTypeEnum.Color3D:
                    if (mDrawerColorMap == null)
                    {
                        mDrawerColorMap = new DrawerColorMap(this.GraphicsDevice, Content, mVertices, mIndices);
                    }
                    mDrawerColorMap.Flat = false;
                    mDrawerAxis.Flat = false;
                    break;

                case GraphicTypeEnum.Color2D:
                    if (mDrawerColorMap == null)
                    {
                        mDrawerColorMap = new DrawerColorMap(this.GraphicsDevice, Content, mVertices, mIndices);
                    }
                    mDrawerColorMap.Flat = true;
                    mDrawerAxis.Flat = true;
                    CenterView();
                    break;


            }
        }




        /// <summary>
        /// 彩色高度图
        /// </summary>
        DrawerColorMap mDrawerColorMap;
        DrawerContourMap mDrawerContour = null;
        DrawerAxis mDrawerAxis = null;


        //插值顶点
        VertexPositionColorNormal[] mVertices;

        //插值顶点索引
        ushort[] mIndices;
        VertexPositionColorNormal[] mEdgeVertices;
        VertexPositionColorNormal[] mDataPointVertices;

        //2D等高线

        //2D等高区域

        //3D等高线
        //3D等高网格

        //物体对象:



        Matrix matrixRotation = Matrix.CreateRotationY(0f);
        //Stopwatch timer;
        List<PointValue> mMappingData;
        bool mCaptureImage = false;
        Image _Image;

        float mMinValue = float.MaxValue;
        float mMaxValue = float.MinValue;

        //float mImageRadius;

        public Image Image
        {
            get { return _Image; }
            set { _Image = value; }
        }


        Matrix mWorldMatrix;
        Matrix mViewMatrix;
        Matrix mProjectionMatrix;

        public XNAMesh3D()
        {

        }

        bool _IsDataSet = false;
        /// <summary>
        /// 
        /// </summary>
        /// <param name="listPoint"></param>
        /// <param name="radius">硅片的半径,单位与采样点的坐标单位相同</param>
        public void SetData(List<PointValue> listPoint, float radius)
        {
            foreach (PointValue point in listPoint)
            {
                if (point.Value < mMinValue) mMinValue = (float)point.Value;
                if (point.Value > mMaxValue) mMaxValue = (float)point.Value;
            }

            mVertices = MeshGenerator.InterportationKrigingVertexPositionColor(listPoint, radius);
            mIndices = MeshGenerator.CreateIndex();

            mDataPointVertices = MeshGenerator.TranslateDataPoints(listPoint, radius);
            //mImageRadius = 25f;
            mMappingData = listPoint;

            //PrimitiveHelper.ClipByCylinder(ref mVertices, ref mIndices, 23.9f, ref mEdgeVertices);//圆柱形剪栽
            PrimitiveHelper.ClipByCylinder(ref mVertices, ref mIndices, MeshGenerator.GridDimension, ref mEdgeVertices);//圆柱形剪栽
            
            mDrawerContour = null;
            mDrawerColorMap = null;
            mDrawerAxis = null;
            PrepareGraphic();
            _IsDataSet = true;
            this.Visible = true;
        }



        protected override void Initialize()
        {
            //timer = Stopwatch.StartNew();
            // Hook the idle event to constantly redraw our animation.
            //Application.Idle += delegate { Invalidate(); };







        }

        /// <summary>
        /// 绘制2D等高线图.
        /// </summary>
        private void DrawContourLine2D()
        {
            //Vector3 cameraPosition = new Vector3(0, 0, 2.5f);
            //float aspect = GraphicsDevice.Viewport.AspectRatio;

            //计算物体的角度
            //Matrix world = Matrix.CreateTranslation(new Vector3(-25f, -25f, 0f)) * Matrix.CreateRotationX(MathHelper.ToRadians(-90)) * matrixRotation;
            mWorldMatrix = Matrix.CreateTranslation(new Vector3(0f, 0f, 0f)) * matrixRotation;
            //Matrix.CreatePerspectiveFieldOfView(1, aspect, 1, 10);

            mDrawerContour.Draw(mWorldMatrix, mViewMatrix, mProjectionMatrix);
        }



        private void Draw3DHeight()
        {



            //// Spin the triangle according to how much time has passed.
            //float time = (float)timer.Elapsed.TotalSeconds;

            //float yaw = time * 0.4f;
            //float pitch = time * 0.4f;
            //float roll = time * 0.4f;

            //Vector3 cameraPosition = new Vector3(0, 0, 2.5f);

            //float aspect = GraphicsDevice.Viewport.AspectRatio;

            //matrixRotation = Matrix.CreateFromYawPitchRoll(yaw, pitch, roll);// CreateRotationY(0.01f);
            //计算物体的角度
            mWorldMatrix = Matrix.CreateTranslation(new Vector3(0f, 0f, 0f)) * matrixRotation;

            //Matrix.CreatePerspectiveFieldOfView(1, aspect, 1, 10);



            // Draw the current primitive.
            Color color = Color.Red;

            //绘制主图
            mDrawerColorMap.Draw(mWorldMatrix, mViewMatrix, mProjectionMatrix);
        }

        private void CaptureImage()
        {
            //mCaptureImage = true;
            if (mCaptureImage)
            {
                mCaptureImage = false;

                int w = GraphicsDevice.PresentationParameters.BackBufferWidth;
                int h = GraphicsDevice.PresentationParameters.BackBufferHeight;
                using (ResolveTexture2D screenshot = new ResolveTexture2D(GraphicsDevice, w, h, 1, SurfaceFormat.Color))
                {
                    // Grab the screenshot
                    GraphicsDevice.ResolveBackBuffer(screenshot);

                    byte[] textureData = new byte[4 * screenshot.Width * screenshot.Height];
                    screenshot.GetData<byte>(textureData);
                    screenshot.Dispose();

                    System.Drawing.Bitmap bmp = new System.Drawing.Bitmap(
                                   screenshot.Width, screenshot.Height,
                                   System.Drawing.Imaging.PixelFormat.Format32bppArgb
                                 );

                    System.Drawing.Imaging.BitmapData bmpData = bmp.LockBits(
                                   new System.Drawing.Rectangle(0, 0, screenshot.Width, screenshot.Height),
                                   System.Drawing.Imaging.ImageLockMode.WriteOnly,
                                   System.Drawing.Imaging.PixelFormat.Format32bppArgb
                                 );

                    IntPtr safePtr = bmpData.Scan0;
                    System.Runtime.InteropServices.Marshal.Copy(textureData, 0, safePtr, textureData.Length);
                    bmp.UnlockBits(bmpData);

                    _Image = bmp;

                    // just some test output
                    //bmp.Save(@"c:\workbench\smile.bmp", System.Drawing.Imaging.ImageFormat.Bmp);

                }

                if (ImageCaptured != null)
                {
                    ImageCaptured(this, new EventArgs());
                }
            }
        }

        protected override void Draw()
        {
            if (!_IsDataSet)
            {
                this.Visible = false;
                return;
            }
            GraphicsDevice.Clear(Color.White);
            mViewMatrix = Matrix.CreateLookAt(new Vector3(0f, 0f, 150f), new Vector3(0.0f, 0.0f, 0.0f), Vector3.Up);
            //mProjectionMatrix = Matrix.CreatePerspectiveFieldOfView(MathHelper.ToRadians(45), (float)GraphicsDevice.Viewport.Width / (float)GraphicsDevice.Viewport.Height, 1.0f, 10000.0f);//透视投影
            mProjectionMatrix = Matrix.CreateOrthographic((float)GraphicsDevice.Viewport.Width / (float)GraphicsDevice.Viewport.Height * (float)MeshGenerator.GridDimension + 5f, (float)MeshGenerator.GridDimension +5f, 1.0f, 10000.0f);//正交投影

            switch (_GrpahicType)
            {
                case GraphicTypeEnum.Color3D:
                case GraphicTypeEnum.Color2D:
                    Draw3DHeight();
                    break;
                case GraphicTypeEnum.ContourLine2D:
                case GraphicTypeEnum.ContourLine3D:
                case GraphicTypeEnum.ContourZone2D:
                case GraphicTypeEnum.ContourZone3D:
                    DrawContourLine2D();
                    break;
            }

            float yaw = 0, roll = 0, pitch = 0;
            TF300.App.GUI.DatabaseUI.XNALibrary.Helpers.MathHelper.ExtractYarRollPithFromMatrix(matrixRotation, ref yaw, ref roll, ref pitch);
            mDrawerAxis.SetYawRollPitchInfo(yaw, roll, pitch);
            mDrawerAxis.Draw(mWorldMatrix, mViewMatrix, mProjectionMatrix);

        }

        private void InitializeComponent()
        {
            this.SuspendLayout();
            this.ResumeLayout(false);

        }
        protected override void OnMouseClick(MouseEventArgs e)
        {
            mCaptureImage = true;
            base.OnMouseClick(e);
        }

        bool mDraging = false;
        float _Yaw = 0f;
        float _Pitch = 0f;
        float _Roll = 0f;
        float _MouseSpeed = 0.01f;

        int _OriginX;
        int _OriginY;

        /// <summary>
        /// 使视图居中
        /// </summary>
        public void CenterView()
        {
            _Yaw = 0f;
            _Pitch = 0f;
            _Roll = 0f;
            matrixRotation = Matrix.CreateFromYawPitchRoll(0f, 0f, 0f);
            lastRationMatrix = Matrix.CreateFromYawPitchRoll(0f, 0f, 0f);
            this.Invalidate();
        }

        /// <summary>
        /// 鼠标按下:设置鼠标捕获并做标记
        /// </summary>
        /// <param name="e"></param>
        protected override void OnMouseDown(MouseEventArgs e)
        {
            mDraging = true;
            this.Capture = true;
            base.OnMouseDown(e);
            _OriginX = e.X;
            _OriginY = e.Y;

        }

        protected override void OnMouseUp(MouseEventArgs e)
        {
            mDraging = false;
            this.Capture = false;
            //lastRationMatrix = lastRationMatrix * Matrix.CreateFromYawPitchRoll(_Yaw, _Pitch, 0f);
            //_Yaw = 0f; _Pitch = 0f;
            base.OnMouseUp(e);
        }
        Matrix lastRationMatrix = Matrix.CreateFromYawPitchRoll(0f, 0f, 0f);
        protected override void OnMouseMove(MouseEventArgs e)
        {
            if (mDraging)
            {

                if (e.Button == MouseButtons.Left && (_LockRoll == false) && (_GrpahicType == GraphicTypeEnum.Color3D || _GrpahicType == GraphicTypeEnum.ContourLine3D || _GrpahicType == GraphicTypeEnum.ContourZone3D))//yaw pitch
                {
                    _Yaw += _MouseSpeed * (e.X - _OriginX);
                    _Pitch += _MouseSpeed * (e.Y - _OriginY);
                }
                else if (e.Button == MouseButtons.Right || _LockRoll || (!(_GrpahicType == GraphicTypeEnum.Color3D || _GrpahicType == GraphicTypeEnum.ContourLine3D || _GrpahicType == GraphicTypeEnum.ContourZone3D)))//Roll
                {
                    _Roll += _MouseSpeed * (e.X - _OriginX);
                }
                this.Invalidate();
                _OriginX = e.X;
                _OriginY = e.Y;
                Matrix RollMatrix = Matrix.CreateFromYawPitchRoll(0f, 0f, _Roll);
                matrixRotation = RollMatrix * lastRationMatrix * Matrix.CreateFromYawPitchRoll(_Yaw, _Pitch, 0f);// CreateRotationY(0.01f);

                lastRationMatrix = lastRationMatrix * Matrix.CreateFromYawPitchRoll(_Yaw, _Pitch, 0f);
                _Yaw = 0f; _Pitch = 0f;


            }

            base.OnMouseMove(e);
        }


        bool _LockRoll = false;
        /// <summary>
        /// 获取或设置鼠标拖动时是否锁定图形绕Y轴旋转
        /// </summary>
        public bool LockRoll
        {
            get { return _LockRoll; }
            set { _LockRoll = value; }
        }



    }






    public enum GraphicTypeEnum
    {
        ContourLine2D,
        ContourLine3D,
        ContourZone2D,
        ContourZone3D,
        Color2D,
        Color3D,

    }


}
