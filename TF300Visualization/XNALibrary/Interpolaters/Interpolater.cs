using System;
using System.Collections.Generic;
using System.Text;
using TF300.App.GUI.DatabaseUI.XNALibrary.Interpolaters;

namespace XNAHelper.Interpolaters
{
    public interface Interpolater
    {
        double GetInterpolatedZ(float xpos, float ypos);
    }

}
