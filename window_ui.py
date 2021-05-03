import pymel.core as pm
import maya.cmds as cmds
from maya import cmds , OpenMaya
import math


# clean channel box
def hideShapes(*args):
    obj = cmds.ls(tr=1, s=1)

    for o in obj:
        cmds.setAttr(o + '.isHistoricallyInteresting', 0)
def hideInputs(*args):
    obj = cmds.ls(dep=0, nt=0, tr=0,s=0,lt=0)

    for o in obj:

        cmds.setAttr(str(o) + '.isHistoricallyInteresting', 0)

# -------------transform objects---------------------------#
# snap objects
def snapObject(*args):
    sel = pm.ls(sl=1)

    target = sel[0]
    ctrlObject = sel[1:]
    for c in ctrlObject:
        prntConst = pm.parentConstraint(target, c, w=1)
        pm.delete(prntConst)
# centerPivot
def centerPivot(target):
    sel = pm.ls(sl=1)
    for s in sel:
        pm.xform(s, cp=1)
def zeroGrp(obj):
    grp_name = str(obj[0]) + '_offset'

    offset_grp = pm.group(em=True, name=grp_name)

    # pos = pm.xform(offset_grp, p = 1)

    prntConst = pm.parentConstraint(obj, offset_grp, mo=0)

    pm.delete(prntConst)
    pm.parent(obj, offset_grp)
    return offset_grp
def snap(obj, target):
    prntConst = pm.parentConstraint(target, obj, mo=0)

    pm.delete(prntConst)
# reset transform
def resetTranslate(*args):
    sel = pm.ls(sl=1)
    for s in sel:
        # pm.connectAttr( 'l_stretch_mult_scaleY.outputX' ,s+'.scaleY', f=1)
        pm.setAttr(s + '.translateX', 0)
        pm.setAttr(s + '.translateY', 0)
        pm.setAttr(s + '.translateZ', 0)
def resetRotation(*args):
    sel = pm.ls(sl=1)
    for s in sel:
        # pm.connectAttr( 'l_stretch_mult_scaleY.outputX' ,s+'.scaleY', f=1)
        pm.setAttr(s + '.rotateX', 0)
        pm.setAttr(s + '.rotateY', 0)
        pm.setAttr(s + '.rotateZ', 0)
def resetScale(*args):
    sel = pm.ls(sl=1)
    for s in sel:
        # pm.connectAttr( 'l_stretch_mult_scaleY.outputX' ,s+'.scaleY', f=1)
        pm.setAttr(s + '.scaleX', 1)
        pm.setAttr(s + '.scaleY', 1)
        pm.setAttr(s + '.scaleZ', 1)

# -------------create heirarchy---------------------------#
# create offsetGrp


def createOffsetGrp(*args):

    target = pm.ls(sl=1)
    ctrl_parent = pm.listRelatives(target, p=1)
    offset_grp = pm.group(target, n=target[0] + 'offset')
    pm.xform(offset_grp, cp=1)
    return offset_grp
def createCtrlLocHierarchy(*args):
    sel = pm.ls(sl=1)

    for s in sel:
        pos = pm.xform(s, q=1, ws=1, t=1)
        rot = pm.xform(s, q=1, ws=1, ro=1)

        loc = pm.spaceLocator()
        grp = pm.group(em=1)
        ctrl = pm.circle(r=1)[0]

        pm.rename(loc, s + "_LOC")
        pm.rename(grp, s + "_grp")
        pm.rename(ctrl, s + "_ctrl")

        pm.xform(loc, ws=1, t=pos)
        pm.xform(loc, ws=1, ro=rot)

        pm.xform(grp, ws=1, t=pos)
        pm.xform(grp, ws=1, ro=rot)

        pm.xform(ctrl, ws=1, t=pos)
        pm.xform(ctrl, ws=1, ro=rot)

        pm.parent(s, loc)

        pm.parent(loc, ctrl)
        pm.parent(ctrl, grp)
def createLocJntHierarchy(*args):
    sel = pm.ls(sl=1)
    print sel
    for s in sel:
        pos = pm.xform(s, q=1, ws=1, t=1)
        rot = pm.xform(s, q=1, ws=1, ro=1)

        loc = pm.spaceLocator()
        grp = pm.group(em=1)

        pm.rename(loc, s + "_LOC")
        pm.rename(grp, s + "_grp")

        pm.xform(loc, ws=1, t=pos)
        pm.xform(loc, ws=1, ro=rot)

        pm.xform(grp, ws=1, t=pos)
        pm.xform(grp, ws=1, ro=rot)

        pm.parent(s, loc)
        pm.parent(loc, grp)

# -------------create joints---------------------------#
'''
# create joints on Crv
def jntsAlongCrv(curve=None, segments=2):
    # creates an oriented joint chain that follows a selected curve
    if segments < 2:
        return
    if not curve:
        curve = cmds.ls(sl=1, o=1)[0]

    # use pointOnCurve(pp=0->2) to find world coordinates on the curve
    jnts = []
    a = 0.0
    # find the max U value. It's equal to the number of spans on the curve
    u = cmds.getAttr(curve + ".spans")
    # these have to be float values!
    i = u / (float(segments) - 1.0)
    cmds.select(cl=1)
    # c will be used to check for computational errors because on some numbers we end at segments - 1
    c = 0
    while a <= u:
        jnt = cmds.joint(p=cmds.pointOnCurve(curve, pr=a))
        pm.parent(jnt, w=1)
        jnts.append(jnt)
        a += i
        c += 1

    # in some cases (segments=10,12,19) we end at one joint before the end. This fixes this problem.
    if c != segments:
        jnts.append(cmds.joint(p=cmds.pointOnCurve(curve, pr=u)))

    return jnts
'''
# non zeo joints
def nonZeroJoints(*args):
    sel = pm.ls(sl=1)

    trans = ['translate', 'rotate', 'scale']

    axis = ['x', 'y', 'z']

    for s in sel:
        rotX = pm.getAttr(s + '.rotateX')
        rotY = pm.getAttr(s + '.rotateY')
        rotZ = pm.getAttr(s + '.rotateZ')

        jntOrntX = pm.getAttr(s + '.jointOrientX')
        jntOrntY = pm.getAttr(s + '.jointOrientY')
        jntOrntZ = pm.getAttr(s + '.jointOrientZ')
        print jntOrntY

        if (rotX != 0):
            jntOrntX = rotX + jntOrntX
            pm.setAttr(s + '.jointOrientX', jntOrntX)

        if (rotY != 0):
            jntOrntY = rotY + jntOrntY
            pm.setAttr(s + '.jointOrientY', jntOrntY)

        if (rotZ != 0):
            jntOrntZ = rotZ + jntOrntZ
            pm.setAttr(s + '.jointOrientZ', jntOrntZ)

        pm.setAttr(s + '.rotateX', 0)
        pm.setAttr(s + '.rotateY', 0)
        pm.setAttr(s + '.rotateZ', 0)
# create joints on the center of selected vertices
def createJointOnCenter(*args):
    # create joint on center
    sel = cmds.ls(sl=1, fl=1)
    xList = []
    yList = []
    zList = []

    for s in sel:
        pos = cmds.xform(s, ws=1, q=1, t=1)
        print pos

        xList.append(pos[0])
        yList.append(pos[1])
        zList.append(pos[2])

    sum_x = sum(xList) / len(xList)
    sum_y = sum(yList) / len(xList)
    sum_z = sum(zList) / len(xList)
    print sum_x

    print sum_y

    print sum_z

    jnt = cmds.joint(p=(sum_x, sum_y, sum_z))
    cmds.parent(jnt, w=1)
# create joints on selection
def createJointOnSelectedVeritces(*args):
    sel = cmds.ls(sl=1, fl=1)
    pm.select(cl=1)
    for s in sel:
        pos = pm.xform(s, ws=1, q=1, t=1)
        jnt = pm.joint(p = pos)
        pm.parent(jnt, w=1)

#zero out joint orientation
def zeroOutJointOrientation(*args):
    sel = pm.ls(sl=1)

    for s in sel:

        jntOrntX = pm.getAttr(s + '.jointOrientX')
        jntOrntY = pm.getAttr(s + '.jointOrientY')
        jntOrntZ = pm.getAttr(s + '.jointOrientZ')

        pm.setAttr(s + '.jointOrientX', 0)
        pm.setAttr(s + '.jointOrientY', 0)
        pm.setAttr(s + '.jointOrientZ', 0)

#create heirarchy of selected joints
def createJointHierarchyOnSelection(*args):
    sel = cmds.ls(sl=1)
    selectionLength = len(sel) - 1
    i = 0
    for i in range(selectionLength):
        pm.parent(sel[i], sel[i + 1])
        i = i + 1

# remove joints
def removeJoints(*args):
    sel = cmds.ls(sl=1)

    for s in sel:
        cmds.removeJoint(s)

# -------------create rig setup---------------------------#
def colorCtr(color=None, ctrl=None, *args):
    if (color == "RED"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(1.0, 0.0, 0.0)

    if (color == "BLUE"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(0.0, 0.0, 1.0)

    if (color == "LIGHTBLUE"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(0.0, 0.156, 1.0)

    if (color == "GREEN"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(0.0, 1.0, 0.0)

    if (color == "YELLOW"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(1.0, 1.0, 0.0)

    if (color == "PINK"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(1.0, 0.0, 1.0)

    if (color == "VOILET"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(0.188, 0.0, 1.0)

    if (color == "ORANGE"):
        ctrl.overrideEnabled.set(1)
        ctrl.overrideRGBColors.set(1)
        ctrl.overrideColorRGB.set(1.0, 0.094, 0.0)
def snap(target, ctrlObject):
    prntConst = pm.parentConstraint(target, ctrlObject, w=1)

    pm.delete(prntConst)
def createOffset(target, name):
    pm.select(clear=1)
    ctrl_parent = pm.listRelatives(target, p=1)
    offset_grp = pm.group(target, n=target[0] + name)
    pm.xform(offset_grp, cp=1)

    return offset_grp
# create jiggle joints
def createJiggleonJoints(*args):
    # pm.parent(obj,offset_grp)

    # get list of the joints
    sel = pm.ls(sl = 1)

    if pm.objExists("tuner_ctrl"):
        tuner_ctrl = "tuner_ctrl"

    else:
        tuner_ctrl = pm.circle(n = "tuner_ctrl")
        pm.addAttr(tuner_ctrl, ln= "Stiffness",at= "double", dv= 1,k=1)

    if pm.objExists("dynamics_grp"):
        dynamics_grp = "dynamics_grp"

    else:
        dynamics_grp = pm.group(n = "dynamics_grp", w=1)
    if pm.objExists("drive_ctrl"):
        drive_ctrl = "drive_ctrl"

    else:
        drive_ctrl= pm.circle(n = "drive_ctrl")


    drive_ctrl = "drive_ctrl"
    #in loop

    for s in sel:

        jiggle_ctrl = pm.circle( n = s+ "_jiggle", nr = [1, 0, 0])

        #create plane
        goal_plane = pm.polyPlane(n = s + "_goal_plane", w=1 , h =1,sx= 1,sy= 1,ax =[ 1, 0, 0],cuv =2 ,ch =1)

        #move it a bit forward
        chld = pm.listRelatives(s, c = 1)

        ref = pm.duplicate(chld)

        plane_move_location = pm.getAttr(ref[0]+".translateX")

        pm.setAttr(ref[0]+".translateX", plane_move_location+1)

        #oreint plane to joint
        snap(goal_plane,ref[0])

        pm.delete(ref)

        #create zero group to zero out the tranlsation
        goal_plane_offset_grp = zeroGrp(goal_plane)

        #create soft body
        soft_plane_particle = pm.soft(goal_plane,c= 0, d= 1, g =0.5)
        soft_plane_particle_shape = soft_plane_particle[0] + "Shape"

        pm.select("copyOf"+ goal_plane[0])
        soft_plane = pm.ls(sl =1)

        #create locator
        loc = pm.spaceLocator(n = s+ "_loc_aim")

        #point on poly constraint locator it on the plane
        pm.pointOnPolyConstraint(soft_plane, loc, offset= [0 ,0, 0],weight=1)

        snap(jiggle_ctrl,goal_plane)
        jiggle_ctrl_grp = zeroGrp(jiggle_ctrl)

        pm.parentConstraint( jiggle_ctrl, goal_plane_offset_grp, mo =0 )

        #create aim constraint with maintain offset for joint to locator
        pm.aimConstraint(loc, s, mo=1,weight= 1,aimVector=[1, 0, 0],upVector=[0, 1, 0],worldUpType= "vector",worldUpVector= [0, 1, 0])

        stiffness_mult = pm.shadingNode("multiplyDivide", asUtility =True)

        pm.connectAttr(str(tuner_ctrl[0])+".Stiffness", str(stiffness_mult)+".input1X", f=1)

        pm.setAttr(str(stiffness_mult)+".input2X", 0.1)

        pm.connectAttr( str(stiffness_mult)+".outputX", str(soft_plane_particle_shape)+".goalWeight[0]", f=1)

        pm.parent(goal_plane_offset_grp,dynamics_grp)
        pm.parent(soft_plane,dynamics_grp)
        pm.parent(loc,dynamics_grp)
        pm.parent(jiggle_ctrl_grp,drive_ctrl)
        pm.parent(sel,drive_ctrl)
# movable pivot on selected control
def createMovablePivot(*args):
    sel = pm.ls(sl=1)
    ctrl = sel[0]
    jnt= sel[1]

    pos = pm.xform(ctrl, ws=1,q=1, t=1)

    pivot = pm.spaceLocator(n="pivot", p= pos)

    ctrl_loc = pm.spaceLocator(n="ctrl_loc", p= pos)

    pm.parent(pivot,ctrl)
    pm.parent(ctrl_loc, ctrl)
    pm.parentConstraint(ctrl_loc, jnt, mo=1)
    pm.setAttr(ctrl_loc+'.visibility',0)

    pm.connectAttr( pivot+'.translateX',ctrl+'.rotatePivotX')
    pm.connectAttr( pivot+'.translateY',ctrl+'.rotatePivotY')
    pm.connectAttr( pivot+'.translateZ',ctrl+'.rotatePivotZ')
# create splineIK on selected vertices

def createFKChain(*args):

    fk_chain = pm.ls(sl=1)

    i = 0
    # sel = pm.ls(selection = True)
    fk_ctrl_list = []
    fk_ctrl_offset_list = []

    for s in fk_chain:
        rotX = pm.getAttr(s + '.rotateX')
        rotY = pm.getAttr(s + '.rotateY')
        rotZ = pm.getAttr(s + '.rotateZ')

        jntOrntX = pm.getAttr(s + '.jointOrientX')
        jntOrntY = pm.getAttr(s + '.jointOrientY')
        jntOrntZ = pm.getAttr(s + '.jointOrientZ')

        if (rotX != 0):
            jntOrntX = rotX + jntOrntX
            pm.setAttr(s + '.jointOrientX', jntOrntX)

        if (rotY != 0):
            jntOrntY = rotY + jntOrntY
            pm.setAttr(s + '.jointOrientY', jntOrntY)

        if (rotZ != 0):
            jntOrntZ = rotZ + jntOrntZ
            pm.setAttr(s + '.jointOrientZ', jntOrntZ)

        pm.setAttr(s + '.rotateX', 0)
        pm.setAttr(s + '.rotateY', 0)
        pm.setAttr(s + '.rotateZ', 0)

        # reorient last joint
    pm.setAttr(fk_chain[-1] + ".jointOrientX", 0)
    pm.setAttr(fk_chain[-1] + ".jointOrientY", 0)
    pm.setAttr(fk_chain[-1] + ".jointOrientZ", 0)

    for i in range(len(fk_chain)):
        fk_ctrl = pm.circle(n=fk_chain[i] + '_fk_ctrl', c=(0, 0, 0), nr=(0, 0, 0), sw=360, r=1, d=3, ut=0, tol=0.01,
                            s=8)
        fk_ctrl_offset = createOffset(fk_ctrl, "offset")
        snap(fk_chain[i], fk_ctrl_offset)

        pm.orientConstraint(fk_ctrl, fk_chain[i], mo=1)

        fk_ctrl_list.append(fk_ctrl)
        fk_ctrl_offset_list.append(fk_ctrl_offset)

        colorCtr("RED", fk_ctrl[0])
        pm.select(clear=1)
        i = i + 1

    for i in range(len(fk_chain) - 1):
        pm.parent(fk_ctrl_offset_list[i + 1], fk_ctrl_list[i][0])


def createIK_chain(guideList, *args):
    # loc
    # guideList= createLimbGuide[0]
    # print(guideList)
    refresh_list = []
    jnt_chain = pm.ls(sl=1)

    root_ctrl = pm.circle(n="root_ctrl", c=(0, 0, 0), nr=(0, 1, 0), sw=360, r=2, d=1, ut=0,
                          tol=0.01, s=8)
    root_ctrl_offset = createOffset(root_ctrl, "offset_grup")
    pm.select(clear=1)

    pos_upper = pm.xform(jnt_chain[0], ws=1, t=1, q=1)
    pos_middle = pm.xform(jnt_chain[1], ws=1, t=1, q=1)
    pos_end = pm.xform(jnt_chain[2], ws=1, t=1, q=1)

    upper_loc = pm.spaceLocator(p=pos_upper)  # upperlimb
    middle_loc = pm.spaceLocator(p=pos_middle) # middlelimb
    end_loc = pm.spaceLocator(p=pos_end) # endlimb

    # control list
    pm.select(clear=1)

    # orient joint chain
    pm.joint(jnt_chain[0], e=1, oj='yzx', secondaryAxisOrient="zup", ch=1, zso=1)
    for s in jnt_chain:
        rotX = pm.getAttr(s + '.rotateX')
        rotY = pm.getAttr(s + '.rotateY')
        rotZ = pm.getAttr(s + '.rotateZ')

        jntOrntX = pm.getAttr(s + '.jointOrientX')
        jntOrntY = pm.getAttr(s + '.jointOrientY')
        jntOrntZ = pm.getAttr(s + '.jointOrientZ')
        print jntOrntY

        if (rotX != 0):
            jntOrntX = rotX + jntOrntX
            pm.setAttr(s + '.jointOrientX', jntOrntX)

        if (rotY != 0):
            jntOrntY = rotY + jntOrntY
            pm.setAttr(s + '.jointOrientY', jntOrntY)

        if (rotZ != 0):
            jntOrntZ = rotZ + jntOrntZ
            pm.setAttr(s + '.jointOrientZ', jntOrntZ)

        pm.setAttr(s + '.rotateX', 0)
        pm.setAttr(s + '.rotateY', 0)
        pm.setAttr(s + '.rotateZ', 0)

    # reorient last joint
    pm.setAttr(jnt_chain[2] + ".jointOrientX", 0)
    pm.setAttr(jnt_chain[2] + ".jointOrientY", 0)
    pm.setAttr(jnt_chain[2] + ".jointOrientZ", 0)

    pm.select(clear=1)
    # ------------------------------------------------------------------------
    # create IK handle
    ik_hdl = pm.ikHandle(sj=jnt_chain[0], ee=jnt_chain[2], n="l_arm_ikHdl")

    #ik_hdl_grp = createOffset(ik_hdl[0],"offset_grup")
    # create controls
    ik_ctrl = pm.circle(n="ik_ctrl")

    pm.rename(ik_ctrl[0], 'l_arm_ik_ctrl')
    snap(jnt_chain[2], ik_ctrl)
    ik_ctrl_grp = createOffset(ik_ctrl,"offset_grup")
    pm.parent(ik_ctrl_grp, root_ctrl[0])

    # constraint ik hdl
    pm.pointConstraint(ik_ctrl, ik_hdl[0], offset=[0, 0, 0], weight=1)
    # -------------------------------------------------------------------------
    # create locators for pole vector position

    loc1 = pm.spaceLocator(n="loc1")
    loc2 = pm.spaceLocator(n="loc2")

    grp1 = pm.group(loc1, n=loc1 + "grp1_offset")
    grp2 = pm.group(loc2, n=loc2 + "grp2_offset")
    # create elbow control
    elbow_ctrl = pm.circle(n="elbow_ctrl")

    pm.rename(elbow_ctrl[0], 'elbow_ctrl')
    elbow_ctrl_grp = createOffset(elbow_ctrl,"offset_grup")
    pm.parent(elbow_ctrl_grp, root_ctrl[0])

    # point constraint grp1 to root and end
    pm.pointConstraint(jnt_chain[0], jnt_chain[2], grp1, offset=[0, 0, 0], weight=1)
    # aim const grp1 to root jnt
    pm.aimConstraint(jnt_chain[0], grp1, offset=[0, 0, 0], weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="scene")
    # point const loc1 to middle joint, only y axis
    pm.pointConstraint(jnt_chain[1], loc1, offset=[0, 0, 0], skip=["x", "z"], weight=1)
    # point const grp2 to middle joint
    pm.pointConstraint(jnt_chain[1], grp2, offset=[0, 0, 0], weight=1)
    # aim constraint grp2 to loc1, choose axis which will move loc1 away
    pm.aimConstraint(loc1, grp2, offset=[0, 0, 0], weight=1, aimVector=[0, 0, 1], upVector=[0, 1, 0],
                     worldUpType="scene")
    pm.move(0, 0, -12, loc2, os=True, rpr=True)
    snap(loc2, elbow_ctrl_grp)

    pm.poleVectorConstraint(elbow_ctrl, ik_hdl[0], w=1)
    # _________________________________________________________________________
    pm.select(clear=1)

    colorCtr("ORANGE", ik_ctrl[0])
    colorCtr("ORANGE", elbow_ctrl[0])

    # delete guides
    pm.delete(grp1)
    pm.delete(grp2)
    pm.delete(upper_loc)
    pm.delete(middle_loc)
    pm.delete(end_loc)

def jntsOnCrv(curve=None, segments=5):
    # creates an oriented joint chain that follows a selected curve
    if segments < 2:
        return
    if not curve:
        curve = cmds.ls(sl=1, o=1)[0]

    # use pointOnCurve(pp=0->2) to find world coordinates on the curve
    jnts = []
    a = 0.0
    # find the max U value. It's equal to the number of spans on the curve
    u = cmds.getAttr(curve + ".spans")
    # these have to be float values!
    i = u / (float(segments) - 1.0)
    cmds.select(cl=1)
    # c will be used to check for computational errors because on some numbers we end at segments - 1
    c = 0
    while a <= u:
        jnt = cmds.joint(p=cmds.pointOnCurve(curve, pr=a))
        pm.parent(jnt, w=1)
        jnts.append(jnt)
        a += i
        c += 1

    # in some cases (segments=10,12,19) we end at one joint before the end. This fixes this problem.
    if c != segments:
        jnts.append(cmds.joint(p=cmds.pointOnCurve(curve, pr=u)))

    return jnts


'''
# centerPivot
def centerPivot(target):
    pm.xform(target, cp=1)
'''

def createSoftMod(*args):
    #based on softmode setup created by rigmarole studio
    sel = pm.ls(sl=1)
    #get shape node
    selShape = pm.listRelatives(sel,s=1)[0]
    #get softmod
    soft = pm.listConnections(selShape,t="softMod")
    #get transform handle
    softHandle = pm.listConnections(soft, t= "transform")
    #soft = pm.listRelatives(softHandle)
    print softHandle

    soft_ctrl = pm.circle(n = "soft_C_001_CTRL", r=2, ch =0)
    colorCtr("BLUE",soft_ctrl[0])
    soft_grup = pm.group(soft_ctrl, n = "soft_C_001_GRUP")

    bulge_ctrl = pm.circle(n = "bulge_C_001_CTRL",r= 1,ch =0)
    colorCtr("Orange",bulge_ctrl[0])
    bulge_grup = pm.group(bulge_ctrl, n = "bulge_C_001_GRUP")

    soft_loc = pm.spaceLocator(n= "soft_C_001_LOCT")

    pm.parent(bulge_grup,soft_ctrl[0])

    prntConst = pm.parentConstraint(sel,soft_grup, mo =0)
    pm.delete(prntConst)

    pm.parentConstraint(soft_ctrl,soft_loc, mo =1)

    pm.addAttr (bulge_ctrl,ln= "bulge",at ="enum",en= "CONTROL:" , k=1)
    pm.addAttr(bulge_ctrl,ln ="Strength",at ="double",dv= 0, k=1)
    pm.addAttr (bulge_ctrl,ln= "Radius" ,at ="double",dv= 0, k=1)

    pm.connectAttr(str(soft_loc)+".translateX",str(soft[0])+".falloffCenterX")
    pm.connectAttr(str(soft_loc)+".translateY",str(soft[0])+".falloffCenterY")
    pm.connectAttr(str(soft_loc)+".translateZ",str(soft[0])+".falloffCenterZ")
    pm.connectAttr(str(bulge_ctrl[0])+".rotateX", str(softHandle[0])+".rotateX")
    pm.connectAttr(str(bulge_ctrl[0])+".rotateY", str(softHandle[0])+".rotateY")
    pm.connectAttr(str(bulge_ctrl[0])+".rotateZ", str(softHandle[0])+".rotateZ")
    pm.connectAttr(str(bulge_ctrl[0])+".Strength", str(soft[0])+".envelope")
    pm.connectAttr(str(bulge_ctrl[0])+".Radius", str(soft[0])+".falloffRadius")

    pm.connectAttr( str(bulge_ctrl[0])+".translateX", str(softHandle[0])+".translateX")
    pm.connectAttr( str(bulge_ctrl[0])+".translateY", str(softHandle[0])+".translateY")
    pm.connectAttr( str(bulge_ctrl[0])+".translateZ", str(softHandle[0])+".translateZ")

def createSplineIK(*args):
    sel = pm.ls(sl=1, fl=1)

    # clear selection to avoid parenting
    pm.select(clear=1)

    if (pm.objExists("rig_C_001_GRUP") == False):
        rig_grup = pm.group(n="rig_C_001_GRUP")

    if (pm.objExists("bodyDeform_C_001_GRUP") == False):
        bodyDef_grup = pm.group(n="bodyDeform_C_001_GRUP")

    rig_grup = "rig_C_001_GRUP"
    bodyDef_grup = "bodyDeform_C_001_GRUP"

    jntList = []
    jntPosList = []
    bodydeform_obj = []

    # create joints on vrtx
    for s in sel:
        # get position of vertex
        pos = pm.xform(s, ws=1, q=1, t=1)

        # create joint
        jnt = pm.joint()
        pm.parent(jnt, w=1)
        pm.xform(jnt, ws=1, t=pos)

        jntPosList.append(pos)
        jntList.append(jnt)

    # create joint chain
    selectionLength = len(jntList) - 1
    i = 0

    for i in range(selectionLength):
        pm.parent(jntList[i], jntList[i + 1])
        i = i + 1

    pm.reroot(jntList[0])

    # orient chain
    pm.joint(jntList[0], e=1, oj="xyz", ch=1)
    pm.joint(jntList[-1], e=1, oj="none")

    crv = pm.curve(p=jntPosList, d=3)
    # create splineIK curve

    pm.select(cl=1)
    segments = 9
    ctrlJnts = jntsOnCrv(str(crv), segments)

    ctrljntPosList = []

    for c in ctrlJnts:
        # get position of vertex
        pos = pm.xform(c, ws=1, q=1, t=1)

        # jnt pos
        ctrljntPosList.append(pos)

    splineIkCrv = pm.curve(p=ctrljntPosList, d=3, n='splineIkCrv')
    pm.delete(crv)

    ikHdl = pm.ikHandle(n='ikh', sj=str(jntList[0]), ee=str(jntList[-1]), sol='ikSplineSolver', ccv=False, pcv=False,
                        c=splineIkCrv)[0]

    pm.select(ctrlJnts, splineIkCrv)
    pm.skinCluster(ctrlJnts, splineIkCrv, tsb=1)

    # create ctrls
    ctrl_list = []
    aim_list = []
    follow_list = []
    offset_grp_list = []
    for c in ctrlJnts:
        ctrl = pm.circle(c=(0, 0, 0), nr=(0, 1, 0), sw=360, r=1, d=3, ut=0, tol=0.01, s=8, ch=0)

        snap(c, ctrl)

        offset_grp = createOffset(ctrl, "_offset")
        follow_grp = createOffset(ctrl, "_follow")
        aim_grp = createOffset(ctrl, "_aim")

        ctrl_list.append(ctrl)
        aim_list.append(aim_grp)
        offset_grp_list.append(offset_grp)
        follow_list.append(follow_grp)

    i = 0

    for i in range(len(ctrlJnts)):
        pm.parentConstraint(ctrl_list[i], ctrlJnts[i], mo=1)
        i += 1

    pm.pointConstraint(ctrl_list[0], ctrl_list[-1], follow_list[4], mo=1)
    pm.pointConstraint(ctrl_list[0], ctrl_list[4], follow_list[2], mo=1)
    pm.pointConstraint(ctrl_list[4], ctrl_list[-1], follow_list[6], mo=1)

    pm.pointConstraint(ctrl_list[0], ctrl_list[2], follow_list[1], mo=1)
    pm.pointConstraint(ctrl_list[2], ctrl_list[4], follow_list[3], mo=1)
    pm.pointConstraint(ctrl_list[4], ctrl_list[6], follow_list[5], mo=1)
    pm.pointConstraint(ctrl_list[6], ctrl_list[8], follow_list[7], mo=1)

    pm.aimConstraint(ctrl_list[4], aim_list[2], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[0][0]))
    pm.aimConstraint(ctrl_list[4], aim_list[6], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[-1][0]))

    pm.aimConstraint(ctrl_list[2], aim_list[1], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[0][0]))
    pm.aimConstraint(ctrl_list[2], aim_list[3], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[4][0]))
    pm.aimConstraint(ctrl_list[6], aim_list[5], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[4][0]))
    pm.aimConstraint(ctrl_list[6], aim_list[7], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[-1][0]))

    root_ctrl = pm.circle(n="root_ctrl", c=(0, 0, 0), nr=(0, 1, 0), sw=360, r=1, d=3, ut=0, tol=0.01, s=8, ch=0)
    root_ctrl_grup = createOffset(root_ctrl, "_offset")

    # create sine deformer
    sineCrv = pm.duplicate(splineIkCrv)
    pm.rename(sineCrv, 'sineCrv')

    pm.select(sineCrv)

    # pm.nonLinear(type= 'sine',lowBound=1, highBound =1,amplitude= 0, wavelength= 2,dropoff= 0, offset= 0)
    sineDef = pm.nonLinear(type='sine')

    pm.addAttr(root_ctrl, ln="sine", at="enum", en="control:", k=1)
    pm.addAttr(root_ctrl, ln="sineSwitch", at="bool", k=1)
    pm.addAttr(root_ctrl, ln="sineAmplitude", at="double", dv=0, k=1)
    pm.addAttr(root_ctrl, ln="sineWaveLength", at="double", dv=2, k=1)
    pm.addAttr(root_ctrl, ln="sineHighBound", at="double", min=0, max=10, dv=1, k=1)
    pm.addAttr(root_ctrl, ln="sineLowBound", at="double", min=-10, max=0, dv=-1, k=1)
    pm.addAttr(root_ctrl, ln="sineDropOff", at="double", min=-1, max=1, dv=0, k=1)

    # create wave deformer
    waveCrv = pm.duplicate(splineIkCrv)
    pm.rename(waveCrv, 'waveCrv')
    pm.select(waveCrv)
    # pm.nonLinear(type= "wave",minRadius =0,maxRadius= 1,amplitude= 0,wavelength =0,dropoff= 0,offset=0)
    waveDef = pm.nonLinear(type="wave")

    pm.addAttr(root_ctrl, ln="wave", at="enum", en="control:", k=1)
    pm.addAttr(root_ctrl, ln="waveSwitch", at="bool", k=1)
    pm.addAttr(root_ctrl, ln="waveAmplitude", at="double", dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveWaveLength", at="double", dv=2, k=1)
    pm.addAttr(root_ctrl, ln="waveOffset", at="double", min=-10, max=10, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveDropOffPosition", at="double", min=-1, max=1, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveMinRadius", at="double", min=0, max=10, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveMaxRadius", at="double", min=0, max=10, dv=0, k=1)

    bshp = pm.blendShape(sineCrv, waveCrv, splineIkCrv)

    pm.connectAttr(str(root_ctrl[0]) + ".sineSwitch", bshp[0] + ".sineCrv")
    pm.connectAttr(str(root_ctrl[0]) + ".waveSwitch", bshp[0] + ".waveCrv")

    pm.connectAttr(str(root_ctrl[0]) + '.sineAmplitude', str(sineDef[0]) + '.amplitude')
    pm.connectAttr(str(root_ctrl[0]) + '.sineHighBound', str(sineDef[0]) + '.highBound')
    pm.connectAttr(str(root_ctrl[0]) + '.sineLowBound', str(sineDef[0]) + '.lowBound')
    pm.connectAttr(str(root_ctrl[0]) + '.sineWaveLength', str(sineDef[0]) + '.wavelength')
    pm.connectAttr(str(root_ctrl[0]) + '.sineSwitch', str(sineDef[0]) + '.offset')

    pm.connectAttr(str(root_ctrl[0]) + '.waveAmplitude', str(waveDef[0]) + '.amplitude')
    pm.connectAttr(str(root_ctrl[0]) + '.waveDropOffPosition', str(waveDef[0]) + '.dropoffPosition')
    pm.connectAttr(str(root_ctrl[0]) + '.waveMaxRadius', str(waveDef[0]) + '.maxRadius')
    pm.connectAttr(str(root_ctrl[0]) + '.waveMinRadius', str(waveDef[0]) + '.minRadius')
    pm.connectAttr(str(root_ctrl[0]) + '.waveOffset', str(waveDef[0]) + '.offset')
    pm.connectAttr(str(root_ctrl[0]) + '.waveWaveLength', str(waveDef[0]) + '.wavelength')

    pm.parent(offset_grp_list, root_ctrl[0])
    bodydeform_obj.append(jntList)
    bodydeform_obj.append(ctrlJnts)
    bodydeform_obj.append(ikHdl)
    bodydeform_obj.append(splineIkCrv)
    bodydeform_obj.append(sineCrv)
    bodydeform_obj.append(waveCrv)
    bodydeform_obj.append(waveDef)
    bodydeform_obj.append(sineDef)

    pm.parent(bodydeform_obj, bodyDef_grup)

    pm.parent(bodyDef_grup, rig_grup)
    pm.parent(root_ctrl_grup, rig_grup)

    # -------------------------------------------------------------------------------
    sel = pm.ls(sl=1, fl=1)

    # clear selectio to avoid parenting
    pm.select(clear=1)

    if (pm.objExists("rig_C_001_GRUP") == False):
        rig_grup = pm.group(n="rig_C_001_GRUP")
        pm.select(clear=1)

    if (pm.objExists("bodyDeform_C_001_GRUP") == False):
        bodyDef_grup = pm.group(n="bodyDeform_C_001_GRUP")
        pm.select(clear=1)

    rig_grup = "rig_C_001_GRUP"
    bodyDef_grup = "bodyDeform_C_001_GRUP"

    jntList = []
    jntPosList = []
    bodydeform_obj = []

    # create joints on vrtx
    for s in sel:
        # get position of vertex
        pos = pm.xform(s, ws=1, q=1, t=1)

        # create joint
        jnt = pm.joint()
        pm.parent(jnt, w=1)
        pm.xform(jnt, ws=1, t=pos)

        jntPosList.append(pos)
        jntList.append(jnt)

    # create joint chain
    selectionLength = len(jntList) - 1
    i = 0

    for i in range(selectionLength):
        pm.parent(jntList[i], jntList[i + 1])
        i = i + 1

    pm.reroot(jntList[0])

    # orient chain
    pm.joint(jntList[0], e=1, oj="xyz", ch=1)
    pm.joint(jntList[-1], e=1, oj="none")

    crv = pm.curve(p=jntPosList, d=3)
    # create splineIK curve

    pm.select(cl=1)
    segments = 9
    ctrlJnts = jntsOnCrv(str(crv), segments)

    ctrljntPosList = []

    for c in ctrlJnts:
        # get position of vertex
        pos = pm.xform(c, ws=1, q=1, t=1)

        # jnt pos
        ctrljntPosList.append(pos)

    splineIkCrv = pm.curve(p=ctrljntPosList, d=3, n='splineIkCrv')
    pm.delete(crv)

    ikHdl = pm.ikHandle(n='ikh', sj=str(jntList[0]), ee=str(jntList[-1]), sol='ikSplineSolver', ccv=False, pcv=False,
                        c=splineIkCrv)[0]

    pm.select(ctrlJnts, splineIkCrv)
    pm.skinCluster(ctrlJnts, splineIkCrv, tsb=1)

    # create ctrls
    ctrl_list = []
    aim_list = []
    follow_list = []
    offset_grp_list = []
    for c in ctrlJnts:
        ctrl = pm.circle(c=(0, 0, 0), nr=(0, 1, 0), sw=360, r=1, d=3, ut=0, tol=0.01, s=8, ch=0)

        snap(c, ctrl)

        offset_grp = createOffset(ctrl, "_offset")
        follow_grp = createOffset(ctrl, "_follow")
        aim_grp = createOffset(ctrl, "_aim")

        ctrl_list.append(ctrl)
        aim_list.append(aim_grp)
        offset_grp_list.append(offset_grp)
        follow_list.append(follow_grp)

    i = 0

    for i in range(len(ctrlJnts)):
        pm.parentConstraint(ctrl_list[i], ctrlJnts[i], mo=1)
        i += 1

    pm.pointConstraint(ctrl_list[0], ctrl_list[-1], follow_list[4], mo=1)
    pm.pointConstraint(ctrl_list[0], ctrl_list[4], follow_list[2], mo=1)
    pm.pointConstraint(ctrl_list[4], ctrl_list[-1], follow_list[6], mo=1)

    pm.pointConstraint(ctrl_list[0], ctrl_list[2], follow_list[1], mo=1)
    pm.pointConstraint(ctrl_list[2], ctrl_list[4], follow_list[3], mo=1)
    pm.pointConstraint(ctrl_list[4], ctrl_list[6], follow_list[5], mo=1)
    pm.pointConstraint(ctrl_list[6], ctrl_list[8], follow_list[7], mo=1)

    pm.aimConstraint(ctrl_list[4], aim_list[2], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[0][0]))
    pm.aimConstraint(ctrl_list[4], aim_list[6], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[-1][0]))

    pm.aimConstraint(ctrl_list[2], aim_list[1], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[0][0]))
    pm.aimConstraint(ctrl_list[2], aim_list[3], mo=1, weight=1, aimVector=[0, -1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[4][0]))
    pm.aimConstraint(ctrl_list[6], aim_list[5], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[4][0]))
    pm.aimConstraint(ctrl_list[6], aim_list[7], mo=1, weight=1, aimVector=[0, 1, 0], upVector=[0, 1, 0],
                     worldUpType="objectrotation", worldUpVector=[0, 1, 0], worldUpObject=str(ctrl_list[-1][0]))

    root_ctrl = pm.circle(n="root_ctrl", c=(0, 0, 0), nr=(0, 1, 0), sw=360, r=1, d=3, ut=0, tol=0.01, s=8, ch=0)
    root_ctrl_grup = createOffset(root_ctrl, "_offset")

    # create sine deformer
    sineCrv = pm.duplicate(splineIkCrv)
    pm.rename(sineCrv, 'sineCrv')

    pm.select(sineCrv)

    # pm.nonLinear(type= 'sine',lowBound=1, highBound =1,amplitude= 0, wavelength= 2,dropoff= 0, offset= 0)
    sineDef = pm.nonLinear(type='sine')

    pm.addAttr(root_ctrl, ln="sine", at="enum", en="control:", k=1)
    pm.addAttr(root_ctrl, ln="sineSwitch", at="bool", k=1)
    pm.addAttr(root_ctrl, ln="sineAmplitude", at="double", dv=0, k=1)
    pm.addAttr(root_ctrl, ln="sineWaveLength", at="double", dv=2, k=1)
    pm.addAttr(root_ctrl, ln="sineHighBound", at="double", min=0, max=10, dv=1, k=1)
    pm.addAttr(root_ctrl, ln="sineLowBound", at="double", min=-10, max=0, dv=-1, k=1)
    pm.addAttr(root_ctrl, ln="sineDropOff", at="double", min=-1, max=1, dv=0, k=1)

    # create wave deformer
    waveCrv = pm.duplicate(splineIkCrv)
    pm.rename(waveCrv, 'waveCrv')
    pm.select(waveCrv)
    # pm.nonLinear(type= "wave",minRadius =0,maxRadius= 1,amplitude= 0,wavelength =0,dropoff= 0,offset=0)
    waveDef = pm.nonLinear(type="wave")

    pm.addAttr(root_ctrl, ln="wave", at="enum", en="control:", k=1)
    pm.addAttr(root_ctrl, ln="waveSwitch", at="bool", k=1)
    pm.addAttr(root_ctrl, ln="waveAmplitude", at="double", dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveWaveLength", at="double", dv=2, k=1)
    pm.addAttr(root_ctrl, ln="waveOffset", at="double", min=-10, max=10, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveDropOffPosition", at="double", min=-1, max=1, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveMinRadius", at="double", min=0, max=10, dv=0, k=1)
    pm.addAttr(root_ctrl, ln="waveMaxRadius", at="double", min=0, max=10, dv=0, k=1)

    bshp = pm.blendShape(sineCrv, waveCrv, splineIkCrv, frontOfChain=True)

    pm.connectAttr(str(root_ctrl[0]) + ".sineSwitch", bshp[0] + ".sineCrv")
    pm.connectAttr(str(root_ctrl[0]) + ".waveSwitch", bshp[0] + ".waveCrv")

    pm.connectAttr(str(root_ctrl[0]) + '.sineAmplitude', str(sineDef[0]) + '.amplitude')
    pm.connectAttr(str(root_ctrl[0]) + '.sineHighBound', str(sineDef[0]) + '.highBound')
    pm.connectAttr(str(root_ctrl[0]) + '.sineLowBound', str(sineDef[0]) + '.lowBound')
    pm.connectAttr(str(root_ctrl[0]) + '.sineWaveLength', str(sineDef[0]) + '.wavelength')
    pm.connectAttr(str(root_ctrl[0]) + '.sineSwitch', str(sineDef[0]) + '.offset')

    pm.connectAttr(str(root_ctrl[0]) + '.waveAmplitude', str(waveDef[0]) + '.amplitude')
    pm.connectAttr(str(root_ctrl[0]) + '.waveDropOffPosition', str(waveDef[0]) + '.dropoffPosition')
    pm.connectAttr(str(root_ctrl[0]) + '.waveMaxRadius', str(waveDef[0]) + '.maxRadius')
    pm.connectAttr(str(root_ctrl[0]) + '.waveMinRadius', str(waveDef[0]) + '.minRadius')
    pm.connectAttr(str(root_ctrl[0]) + '.waveOffset', str(waveDef[0]) + '.offset')
    pm.connectAttr(str(root_ctrl[0]) + '.waveWaveLength', str(waveDef[0]) + '.wavelength')

    pm.parent(offset_grp_list, root_ctrl[0])
    bodydeform_obj.append(jntList[0])
    bodydeform_obj.append(ctrlJnts)
    bodydeform_obj.append(ikHdl)
    bodydeform_obj.append(splineIkCrv)
    bodydeform_obj.append(sineCrv)
    bodydeform_obj.append(waveCrv)
    bodydeform_obj.append(waveDef)
    bodydeform_obj.append(sineDef)

    for o in bodydeform_obj:
        pm.parent(o, bodyDef_grup)

    pm.parent(bodyDef_grup, rig_grup)
    pm.parent(root_ctrl_grup, rig_grup)

# Launch the UI
'''This function will house all the commands to build the UI'''

win_name = 'BasicUI'
win_title = 'Basic Window'
win_width = 225
win_height = 75
# If the window already exists, delete it before creating a new one

if pm.window(win_name, exists=True):
    pm.deleteUI(win_name)

# Create the window and its contents
with pm.window(win_name, title=win_title) as win:
    # Declare Controls Here
    with pm.frameLayout(cll=1, label="Transform objects", bgc=[0.3, 0, 0]):
        # -----------------------------------------------------------------------

        with pm.columnLayout(nch=5, adjustableColumn=True, bgc=[0.3, 0, 0]):
            pm.separator(style='in')
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Snap Objects', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select source and "
                                                                                              "target",
                          command=snapObject)
                pm.button(label='Center Pivot', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects",
                          command=centerPivot)
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Reset Translate', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects to "
                                                                                                 "reset their "
                                                                                                 "translation",
                          command=resetTranslate)
                pm.button(label='Reset Rotation', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects to "
                                                                                                "reset their rotation",
                          command=resetRotation)
                pm.button(label='Reset Scale', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects to reset "
                                                                                             "their scale",
                          command=resetScale)
    with pm.frameLayout(cll=1, label="Channel Box", bgc=[0, 0, 0]):
        # -----------------------------------------------------------------------

        with pm.columnLayout(nch=5, adjustableColumn=True, bgc=[0, 0, 0]):
            pm.separator(style='in')
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Hide Shapes', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Cleanes channel box by "
                                                                                             "hiding shaped of "
                                                                                             "transforms",
                          command= hideShapes)
                pm.button(label='Hide Inputs', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Celanes channel box by "
                                                                                             "hidiing inputs of "
                                                                                             "transform",
                          command=hideInputs)

    with pm.frameLayout(cll=1, label="Create Hierarchy", bgc=[0.2, 0.3, 0.2]):
        with pm.columnLayout(nch=3, adjustableColumn=True, bgc=[0.2, 0.3, 0.2]):
            pm.separator(style='in')
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Offset Grup', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects",
                          c=createOffsetGrp)
                pm.button(label='CtrlLoc', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects",
                          c=createCtrlLocHierarchy)
                pm.button(label='LocJnt', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select objects",
                          c=createLocJntHierarchy)

    with pm.frameLayout(cll=1, label="Joints Operation", bgc=[0, 0, 0.2]):
        with pm.columnLayout(nch=5, adjustableColumn=True, bgc=[0, 0, 0.2]):
            pm.separator(style='in')
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Vertices', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select vertices",
                          c=createJointOnSelectedVeritces)
                pm.button(label='Center of Selection', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Select "
                             "vertices to "
                             "create join "
                             "on center",c=createJointOnCenter)
                pm.button(label='Remove Joints', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Select joints to "
                               "remove them from "
                               "chain",
                          c=removeJoints)
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left'), h=win_height / 2):
                pm.button(label='NonZero', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Cleans rotation value of joints by adding them to joint orientation", c=nonZeroJoints)
                pm.button(label='ZeroJointOrient', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Sets the joint "
                             "orientation to "
                             "zero", c=zeroOutJointOrientation)
                pm.button(label='Joint Heirarchy', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Create heirarchy in selected order", c= createJointHierarchyOnSelection)

    with pm.frameLayout(cll=1, label="Create rig setups", bgc=[0, 0.2, 0.2]):
        with pm.columnLayout(nch=6, adjustableColumn=True, bgc=[0, 0.2, 0.2]):
            pm.separator(style='in')
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='IK chain', ekf=1, w=win_width / 2, h=win_height / 2,
                          sbm="Select 3 joints in heirarchy",c =createIK_chain )
                pm.button(label='FK chain', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select joints in "
                                                                                          "heirarchy", c= createFKChain)
                pm.button(label='Bulge', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Create bulge control on "
                                                                                       "geometry with softdeformer",
                          c=createSoftMod)
            with pm.rowLayout(nc=3, adjustableColumn=True, cal=(1, 'left')):
                pm.button(label='Movable Pivot', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select ctrl and joint "
                                                                                               "to be constrained", c=createMovablePivot)
                pm.button(label='Jiggle Joints', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select joints",
                          c= createJiggleonJoints)
                pm.button(label='Motion Path Joint', ekf=1, w=win_width / 2, h=win_height / 2, sbm="Select curve and "
                                                                                                   "joints")

# Show the window and edit its size

pm.showWindow(win_name)
pm.window(win_name, edit=True, widthHeight=(win_width, win_height), rtf=1)

