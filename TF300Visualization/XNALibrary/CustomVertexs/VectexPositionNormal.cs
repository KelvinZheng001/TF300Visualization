using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework;

namespace XNAHelper.CustomVertexs
{

    /// <summary>
    /// Custom vertex type for vertices that have just a
    /// position and a normal, without any texture coordinates.
    /// </summary>
    public struct VertexPositionNormal
    {
        public Vector3 Position;
        public Vector3 Normal;


        /// <summary>
        /// Constructor.
        /// </summary>
        public VertexPositionNormal(Vector3 position, Vector3 normal)
        {
            Position = position;
            Normal = normal;
        }


        /// <summary>
        /// Vertex format information, used to create a VertexDeclaration.
        /// </summary>
        public static readonly VertexElement[] VertexElements =
        {
            new VertexElement(0, 0, VertexElementFormat.Vector3,
                                    VertexElementMethod.Default,
                                    VertexElementUsage.Position, 0),

            new VertexElement(0, 12, VertexElementFormat.Vector3,
                                     VertexElementMethod.Default,
                                     VertexElementUsage.Normal, 0),
        };


        /// <summary>
        /// Size of this vertex type.
        /// </summary>
        public const int SizeInBytes = 24;
    }
}
