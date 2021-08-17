import sys
import cgi, os
import math
import cgitb; cgitb.enable()

import ezdxf
import numpy as np
import copy
import argparse
import PySimpleGUI as sg


from flask import Flask, render_template, request, redirect, request
from werkzeug.utils import secure_filename



app = Flask(__name__)

path=os.getcwd()
UPLOAD_FOLDER=os.path.join(path, 'student')
CORRECT_FOLDER=os.path.join(path, 'correct')
if not os.path.isdir(CORRECT_FOLDER):
    os.mkdir(CORRECT_FOLDER)
app.config['CORRECT_FOLDER']=CORRECT_FOLDER

if not os.path.isdir(UPLOAD_FOLDER):
    os.mkdir(UPLOAD_FOLDER)
app.config['UPLOAD_FOLDER']=UPLOAD_FOLDER


# CHANGE THESE TO RUN ON YOUR SYSTEM
# Locarion of student files to grade
# Location of the correct file


@app.route('/')
def index():
    return render_template('index.html')
@app.route('/my-link/', methods=['POST'])
def main():
    rcorrect_files=[]
    rstudent_files=[]
    raw_student_files =request.files.getlist("student_files")
    raw_correct_files=request.files.getlist("correct_files")
    for file in raw_student_files:
        if file:
            filename=secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
            rstudent_files.append(filename)
    for file in raw_correct_files:
        if file:
            filename=secure_filename(file.filename)
            file.save(os.path.join(app.config['CORRECT_FOLDER'],filename))
            rcorrect_files.append(filename)

    thresh=int(request.form["thresh"])
    extra_ent_penalty = int(request.form["extra_ent_penalty"])
    hatch_error_penalty = int(request.form["hatch_error_penalty"])
    color_error_penalty = int(request.form["color_error_penalty"])
    lw_error_penalty = int(request.form["lw_error_penalty"])
    verbose=int(request.form["verbose"])
    scaling_error=int(request.form["scaling_error"])
    rot_error=int(request.form["rot_error"])
    print(len(rstudent_files))
    correct_files = { filename:ezdxf.readfile(CORRECT_FOLDER + "/" + filename) for filename in rcorrect_files if filename.__contains__(".dxf")}
    student_files = { filename:ezdxf.readfile(UPLOAD_FOLDER + "/" + filename) for filename in rstudent_files if filename.__contains__(".dxf")}
    print(len(student_files))

    
    original_stdout = sys.stdout
    #f = open("GradeCADOutput.txt","w")
    with open ('GradeCADOutput.txt', 'w') as f:
        sys.stdout = f
        for file in student_files:
            for cfile in correct_files:
                print(file, ": ",(grade(student_files[file], correct_files[cfile], verbose, thresh,extra_ent_penalty, hatch_error_penalty, color_error_penalty, lw_error_penalty, scaling_error, rot_error)))
        
    return("Grading succesful. Results have been output in 'GradeCADOutput.txt' file")

# gets the area of a polygon defined by sets of xy coordinates
# xy[:,0] are the x values, and xy[:,1] are the y values
def poly_area(xy):
    l = len(xy)
    s = 0.0
    for i in range(l):
        j = (i+1)%l
        s += 0.5*(xy[i][0]*xy[j][1] - xy[j][0]*xy[i][1])
    return s

# gets the area of the edge paths which define each hatching region
# algorithms used are the same as those used in libreCAD based on taking a line integral to find the area
def get_area(path):
    area = 0
#    area of a polygon
    if path["type"] == "PolyLinePath":
        area += poly_area(path["vertices"])
#    area of a shape defined by lines, arcs, and/or ellipse arcs
    if path["type"] == "EdgePath":
        for edge in path["Edges"]:
            edge = path["Edges"][edge]
            if edge["type"] == "LineEdge":
                area += 0.5*(edge["start"][0]*edge["end"][1] - edge["end"][0]*edge["start"][1])
            if edge["type"] == "ArcEdge":
                a0 = edge["start angle"]*np.pi/180
                a1 = edge["end angle"]*np.pi/180
                r2 = (edge["radius"]**2)/4
                fStart= edge["center"][0]*edge["radius"]*np.sin(a0)+r2*np.sin(2*a0)
                fEnd= edge["center"][0]*edge["radius"]*np.sin(a1)+r2*np.sin(2*a1)
                if edge["is CCW?"] == True:
                    area += fStart - fEnd + 2*r2* abs(a1-a0)
                else:
                    area += fEnd - fStart + 2*r2* abs(a1-a0)
            if edge["type"] == "EllipseEdge":
                majVec = [edge["major axis"][0], edge["major axis"][1]]
                minVec = [edge["minor axis"][0], edge["minor axis"][1]]
                a0 = edge["start angle"]*np.pi/180
                a1 = edge["end angle"]*np.pi/180
                aE = np.arctan(majVec[1]/majVec[0])
                major = np.sqrt( majVec[0]**2 + majVec[1]**2 )
                minor = np.sqrt( minVec[0]**2 + minVec[1]**2 )
                ab = major*minor
                r2 = major**2 + minor**2
                c_x = edge["center"][0]
                start_x = major*minor/np.sqrt(minor**2 + major**2*(np.tan(a0)**2))
                start_x = abs(start_x) if abs(a0) < np.pi/2 else -abs(start_x)
                start_y = np.sqrt(1-(start_x/major)**2)*minor
                start_y = np.sqrt(start_x**2 + start_y**2)*np.sin(aE)+edge["center"][1]
                end_x = major*minor/np.sqrt(minor**2 + major**2*(np.tan(a1)**2))
                end_x = abs(end_x) if abs(a1) < np.pi/2 else -abs(end_x)
                end_y = np.sqrt(1-(end_x/major)**2)*minor
                end_y = np.sqrt(end_x**2 + end_y**2)*np.sin(aE)+edge["center"][1]
                fStart = c_x*start_y+0.25*r2*np.sin(2*aE)*np.cos(a0)*np.cos(a0)-0.25*ab*(2*np.sin(aE)*np.sin(aE)*np.sin(2*a0)-np.sin(2.*a0))
                fEnd = c_x*end_y+0.25*r2*np.sin(2*aE)*np.cos(a1)*np.cos(a1)-0.25*ab*(2*np.sin(aE)*np.sin(aE)*np.sin(2*a1)-np.sin(2*a1))
                if edge["is CCW?"] == True:
                    area += fStart - fEnd + 0.5*ab*(a1-a0)
                else:
                    area += fEnd - fStart + 0.5*ab*(a1-a0)
    return round(area, 5)

# gets the areas of hatch objects defined either by PolyLine Paths or edges which can be lines, arcs, or ellipse arcs
# hatches is the hatch portion of the entities dictionary made by get_ents()
def get_hatch_area(hatches):
    area = 0
    for hatch in hatches:
        for path in hatches[hatch]:
            area += get_area(hatches[hatch][path])
    return round(area, 5)
#caant area polyline hatches
# extract the information about each entity we care about from the dxf file and creates an entity dictionary
def get_ents(dxf):
    header_var_count = len(dxf.header)
    layer_count = len(dxf.layers)
    block_definition_count = len(dxf.blocks)
    entity_count = len(dxf.entities)

    ents = {}
    LWP = [entity for entity in dxf.modelspace() if entity.dxftype() == 'LWPOLYLINE']
    if (LWP != []):
        for i in LWP:
            i.explode()
    Line = [entity for entity in dxf.modelspace() if entity.dxftype() == 'LINE']
    ents["Lines"] = {}
    if (Line != []):
        count = 0
        for i in Line:
            line = {}
            line["layer"] = i.dxf.layer
            line["start"] = [x for x in i.dxf.start]
            line["end"] = [x for x in i.dxf.end]
            line["lineweight"] = i.dxf.lineweight
            line["color"] = i.dxf.color
            ents["Lines"]["Line "+str(count)] = line
            count += 1
    Arc = [entity for entity in dxf.modelspace() if entity.dxftype() == 'ARC']
    ents["Arcs"] = {}
    if (Arc != []):
        count = 0
        for i in Arc:
            arc = {}
            arc["layer"] = i.dxf.layer
            arc["center"] = [x for x in i.dxf.center]
            arc["radius"] = i.dxf.radius
            arc["start angle"] = i.dxf.start_angle
            if (arc["start angle"]>360):
                arc["start angle"]-=360
            if (arc["start angle"]<0):
                arc["start angle"]+=360
            arc["end angle"] = i.dxf.end_angle
            if (arc["end angle"]>360):
                arc["end angle"]-=(360)
            if (arc["end angle"]<0):
                arc["end angle"]+=(360)
            arc["lineweight"] = i.dxf.lineweight
            arc["color"] = i.dxf.color
            ents["Arcs"]["Arc "+str(count)] = arc
            count += 1
    Circ = [entity for entity in dxf.modelspace() if entity.dxftype() == 'CIRCLE']
    ents["Circles"] = {}
    if (Circ != []):
        count = 0
        for i in Circ:
            circ = {}
            circ["layer"] = i.dxf.layer
            circ["center"] = [x for x in i.dxf.center]
            circ["radius"] = i.dxf.radius
            circ["lineweight"] = i.dxf.lineweight
            circ["color"] = i.dxf.color
            ents["Circles"]["Circle "+str(count)] = circ
            count += 1
    Ell = [entity for entity in dxf.modelspace() if entity.dxftype() == 'ELLIPSE']
    ents["Ellipses"] = {}
    if (Ell != []):
        count = 0
        for i in Ell:
            ell = {}
            ell["layer"] = i.dxf.layer
            ell["center"] = [x for x in i.dxf.center]
            ell["major axis"] = [x for x in i.dxf.major_axis]
            ell["minor axis"] = [x for x in i.minor_axis]
            if (ell["major axis"][1]<0):
                ell["major axis"] = [-x for x in i.dxf.major_axis]
            if (ell["minor axis"][1]<0):
                ell["minor axis"] = [-x for x in i.minor_axis]
            ents["Ellipses"]["Ellipse "+str(count)] = ell
            ell["lineweight"] = i.dxf.lineweight
            ell["color"] = i.dxf.color
            count += 1
    Hatch = [entity for entity in dxf.modelspace() if entity.dxftype() == 'HATCH']
    ents["Hatches"] = {}
    if (Hatch != []):
        count = 0
        for i in Hatch:
            hats = {}
            paths = i.paths
            pattern = i.pattern
            paths.polyline_to_edge_path()
            count_path = 0
            for path in paths.paths:
                hats["Path "+str(count_path)] = {}
                if (type(path) == ezdxf.entities.hatch.PolylinePath):
                    hats["Path "+str(count_path)]["type"] = "PolyLinePath"
                    hats["Path " +str(count_path)]["Pattern"] = str(i.dxf.pattern_name)
                    hats["Path "+str(count_path)]["Color"] =str(i.dxf.color)
                    hats["Path "+str(count_path)]["flag"] = path.path_type_flags
                    hats["Path "+str(count_path)]["vertices"] = path.vertices
                else:
                    count_edge = 0
                    hats["Path "+str(count_path)]["type"] = "EdgePath"
                    hats["Path " +str(count_path)]["Pattern"] = str(i.pattern)
                    hats["Path "+str(count_path)]["Color"] =str(i.dxf.color)
                    hats["Path "+str(count_path)]["flag"] = path.path_type_flags
                    edges = {}
                    for edge in path.edges:
                        edges["Edge "+str(count_edge)] = {}
                        if (type(edge) == ezdxf.entities.hatch.LineEdge):
                            edges["Edge "+str(count_edge)]["type"] = "LineEdge"
                            edges["Edge "+str(count_edge)]["start"] = [x for x in edge.start]
                            edges["Edge "+str(count_edge)]["end"] = [x for x in edge.end]
                        if (type(edge) == ezdxf.entities.hatch.ArcEdge):
                            edges["Edge "+str(count_edge)]["type"] = "ArcEdge"
                            edges["Edge "+str(count_edge)]["center"] = [x for x in edge.center]
                            edges["Edge "+str(count_edge)]["radius"] = edge.radius
                            edges["Edge "+str(count_edge)]["start angle"] = edge.start_angle
                            edges["Edge "+str(count_edge)]["end angle"] = edge.end_angle
                            edges["Edge "+str(count_edge)]["is CCW?"] = edge.ccw
                        if (type(edge) == ezdxf.entities.hatch.EllipseEdge):
                            edges["Edge "+str(count_edge)]["type"] = "EllipseEdge"
                            edges["Edge "+str(count_edge)]["center"] = [x for x in edge.center]
                            edges["Edge "+str(count_edge)]["major axis"] = edge.major_axis
                            edges["Edge "+str(count_edge)]["minor axis"] = edge.ratio * edge.major_axis
                            edges["Edge "+str(count_edge)]["start angle"] = edge.start_angle
                            edges["Edge "+str(count_edge)]["end angle"] = edge.end_angle
                            edges["Edge "+str(count_edge)]["is CCW?"] = edge.ccw
                        count_edge += 1
                    hats["Path "+str(count_path)]["Edges"] = edges
                
                    hats["Path "+str(count_path)]["area"] = get_area(hats["Path "+str(count_path)])           
                    count_path += 1
            ents["Hatches"]["Hatch "+str(count)] = hats
            count += 1
    ents["Hatching Area"] = get_hatch_area(ents["Hatches"])
    Lay = [layer for layer in dxf.layers]
    ents["Layers"] = {}
    if (Lay != []):
        count = 0
        for i in Lay:
            layer = {}
            layer["lineweight"] = i.dxf.lineweight
            layer["color"] = i.color
            ents["Layers"][i.dxf.name] = layer
            count += 1
    return ents        

def distl(lin):
    return lin["distance"]
def dista(arc):
    return arc["radius"]
def find_factor(entities_student, entities_correct,thresh):
    factors=[]
    correct_lines = []
    line_factor = 0
    arc_factor = 0
    circ_factor = 0
    ell_factor = 0
    x=0
    for lin_c in entities_correct["Lines"]:
        lin_c = entities_correct["Lines"][lin_c]
        distc = (math.dist(lin_c["start"], lin_c["end"]))
        lin_c["distance"] = (distc)
        correct_lines.append(lin_c)
    correct_lines.sort(key = distl)
    student_lines = []
    for lin_s in entities_student["Lines"]:
        lin_s = entities_student["Lines"][lin_s]
        dists = (math.dist(lin_s["start"], lin_s["end"]))
        lin_s["distance"] = dists
        student_lines.append(lin_s)
    student_lines.sort(key = distl)
    while x < len(student_lines) and x < len(correct_lines):
        line_factor = line_factor + (correct_lines[x]['distance']/student_lines[x]['distance'])
        x=x+1
    if x>0:
        line_factor = line_factor/x
    factors.append(line_factor)
    x=0
    correct_arcs = []
    for arc_c in entities_correct["Arcs"]:
        arc_c = entities_correct["Arcs"][arc_c]
        arc_c["radius"] = (arc_c["radius"])
        correct_arcs.append(arc_c)
    correct_arcs.sort(key = dista)
    student_arcs = []
    for arc_s in entities_student["Arcs"]:
        arc_s = entities_student["Arcs"][arc_s]
        arc_s["radius"] = (arc_s["radius"])
        student_arcs.append(arc_s)
    student_arcs.sort(key = dista)
    while x < len(student_arcs) and x < len(correct_arcs):
        arc_factor = arc_factor + (correct_arcs[x]['radius']/student_arcs[x]['radius'])
        x=x+1
    if x>0:
        arc_factor = arc_factor/x
    factors.append(arc_factor)
    x=0
    correct_circs=[]
    for cir_c in entities_correct["Circles"]:
        cir_c = entities_correct["Circles"][cir_c]
        cir_c["radius"] = round(cir_c["radius"], thresh)
        correct_circs.append(cir_c)
    correct_circs.sort(key = dista)
    student_circs=[]
    for cir_s in entities_student["Circles"]:
        cir_s = entities_student["Circles"][cir_s]
        cir_s["radius"] =(cir_s["radius"])
        student_circs.append(cir_s)
    student_circs.sort(key = dista)
    while x < len(student_circs) and x < len(correct_circs):
        circ_factor = circ_factor + (correct_circs[x]['radius']/student_circs[x]['radius'])
        x=x+1
    if x>0:
        circ_factor = circ_factor/x
    factors.append(circ_factor)
    x=0
    correct_ells= []
    for ell_c in entities_correct["Ellipses"]:
        ell_c = entities_correct["Ellipses"][ell_c]
        distc = (math.dist([0,0,0], ell_c["major axis"]))
        ell_c["radius"] = distc
        correct_ells.append(ell_c)
    correct_ells.sort(key = dista)
    student_ells=[]
    for ell_s in entities_student["Ellipses"]:
        ell_s = entities_student["Ellipses"][ell_s]
        ell_s["radius"] = (math.dist([0,0,0], ell_s["major axis"]))
        student_ells.append(ell_s)
    student_ells.sort(key = dista)
    while x < len(student_ells) and x < len(correct_ells):
        ell_factor = ell_factor + (correct_ells[x]['radius']/student_ells[x]['radius'])
        x=x+1
    if x>0:
        ell_factor =ell_factor/x
    factors.append(ell_factor)
    x=0
    return factors
   
def scale(entity_dict, factor):
    ents = copy.deepcopy(entity_dict)
    for lin in ents["Lines"]:
        lin = ents["Lines"][lin]
        lin["start"][0], lin["start"][1] = factor * lin["start"][0], factor * lin["start"][1]
        lin["end"][0], lin["end"][1] = factor * lin["end"][0], factor * lin["end"][1]
    for arc in ents["Arcs"]:
        arc = ents["Arcs"][arc]
        arc["radius"] = factor * arc["radius"]
        arc["center"][0], arc["center"][1] = factor * arc["center"][0], factor * arc["center"][1]
    for cir in ents["Circles"]:
        cir = ents["Circles"][cir]
        cir["center"][0], cir["center"][1] = factor * cir["center"][0], factor * cir["center"][1]
        cir["radius"] = factor * cir["radius"]
    for ell in ents["Ellipses"]:
        ell = ents["Ellipses"][ell]
        ell["center"][0], ell["center"][1] = factor * ell["center"][0], factor * ell["center"][1]
        ell["major axis"][0], ell["major axis"][1] = factor * ell["major axis"][0], factor * ell["major axis"][1]
        ell["minor axis"][0], ell["minor axis"][1] = factor * ell["minor axis"][0], factor * ell["minor axis"][1]
    
    return ents

# The function that takes two entity dictionaries and compares them to produce a grade
def grade_ents(entities_student, entities_correct, factor, verbose, thresh, extra_ent_penalty, hatch_error_penalty, color_error_penalty, lw_error_penalty):
#    the threshold for how similar two values must be to be considered the same
    ent_count = 0
    missing_ents = 0
    ln_count=0
    correct = 0
    extra_ents = 0
    lin_flags=0
    arc_flags=0
    circ_flags=0
    ell_flags=0
    rotation= False
    rotations = 0
    

#    list to ensure each line is counted only once and also to check for extra lines
    
    student_lines = []
    for x in range(0,len(entities_correct["Lines"])):
        lin_c = ((entities_correct["Lines"]["Line " + str(x)]))
        ent_count += 1
        ln_count+=1
        lin_flag = False
#        check if each line from the correct file is found in the student file
        for y in range (0, len(entities_student["Lines"])):
            lin_s = ((entities_student["Lines"]["Line " + str(y)]))
            if (lin_s not in student_lines):
                point1s = [lin_s["start"][0], lin_s["start"][1]]
                point2s= [lin_s["end"][0], lin_s["end"][1]]
                point1c = [lin_c["start"][0], lin_c["start"][1]]
                point2c = [lin_c["end"][0], lin_c["end"][1]]
                dist1 = round(math.dist(point1s, point2s),thresh)
                dist2 = round(math.dist(point1c, point2c),thresh)
                case1 = dist1== dist2
                if x < len((entities_correct["Lines"]))-1 and y < len((entities_student["Lines"]))-1:
                    match1s = [lin_s["end"][0], lin_s["end"][1]]
                    match2s = [entities_student["Lines"]["Line " + str(y+1)]["start"][0], entities_student["Lines"]["Line " + str(y+1)]["start"][1]]
                    match1c = [lin_c["end"][0], lin_c["end"][1]]
                    match2c = [entities_correct["Lines"]["Line " + str(x+1)]["start"][0], entities_correct["Lines"]["Line " + str(x+1)]["start"][1]]
                    if (round(math.dist(match1s, match2s), thresh)) == (round(math.dist(match1c, match2c), thresh)):
                        case2 = True
                    else:
                        case2 = False
                else:
                    case2 = True
                case3a = entities_student["Layers"][lin_s["layer"]]["color"] == entities_correct["Layers"][lin_c["layer"]]["color"]
                case4b = entities_student["Layers"][lin_s["layer"]]["lineweight"] == entities_correct["Layers"][lin_c["layer"]]["lineweight"]
                case3b = lin_s["color"] == lin_c["color"]
                case4a = lin_s["lineweight"] == lin_c["lineweight"]
                if case1  and case2 and (lin_s not in student_lines):
                    if ((point2s[0]-point1s[0])) != 0 and ((point2c[0]-point1c[0])) != 0 and round(((point2s[1] - point1s[1])/(point2s[0]-point1s[0])), thresh) != round(((point2c[1] - point1c[1])/(point2c[0]-point1c[0])), thresh):
                        rotations = rotations + 1
                    correct+=1
                    student_lines.append(lin_s)
                    lin_flag = True
                    if not case3a or not case3b:
                        correct - ((color_error_penalty/100) * 1)
                        color_error = True
                    if not case4a or not case4b:
                        correct - ((lw_error_penalty/100) * 1)
                        lw_error = True
                    break
        if not lin_flag:
            lin_flags = lin_flags+1
#        if the line is missing and verbose is 1, print which line is missing
#        if verbose is 2, then eveything will be printed (mostly for debugging)
        if not lin_flag:
            missing_ents = missing_ents +1
        if verbose!=0:
            if not lin_flag:
                print("missing Line: ",lin_c)
            else:
                if verbose==2:
                    print("found Line: ",lin_c)
#    Now go through each line in the student file and see if it was added to the list of correct lines
#    if there is a line that was not correct, it will be considered extra and points will be deducted
    for lin_s in entities_student["Lines"]:
        lin_s = entities_student["Lines"][lin_s]
        if lin_s not in student_lines:
            extra_ents = extra_ents + 1
            if verbose!=0:
                print("extra Line: ",lin_s)
    while extra_ents != 0 and lin_flags !=0:
        extra_ents = extra_ents -1
        lin_flags = lin_flags -1

#    list to ensure each arc is counted only once and also to check for extra arcs
    student_arcs = []
    for x in range(0, len(entities_correct["Arcs"])):
        arc_c = entities_correct["Arcs"]["Arc " +str(x)]
        ent_count += 1
        arc_flag = False
#        check if each arc from the correct file is also found in the student file
        for y in range(0, len(entities_student["Arcs"])):
            arc_s = entities_student["Arcs"]["Arc " + str(y)]
            if arc_s not in student_arcs:
                case1 = (round(arc_c["end angle"] - arc_c["start angle"], thresh) == round(arc_s["end angle"] - arc_s["start angle"], thresh))
                case2 = round(abs(arc_s["radius"]), thresh) == round(abs(arc_s["radius"]), thresh)
                if len(entities_correct["Lines"]) >0 and len(student_lines) > 0:
                    pointAc = [arc_c["center"][0], arc_c["center"][1]]
                    pointlc = [entities_correct["Lines"]["Line " + str(0)]["end"][0], entities_correct["Lines"]["Line " + str(0)]["end"][1]]
                    pointAs = [arc_s["center"][0], arc_s["center"][1]]
                    pointls = [student_lines[0]["end"][0], student_lines[0]["end"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case3= True
                    else:
                        case3= False
                else:
                    case3 = True            
                case4a = entities_student["Layers"][arc_s["layer"]]["color"] == entities_correct["Layers"][arc_c["layer"]]["color"]
                case5b = entities_student["Layers"][arc_s["layer"]]["lineweight"] == entities_correct["Layers"][arc_c["layer"]]["lineweight"]
                case4b = arc_s["color"] == arc_c["color"]
                case5a = arc_s["lineweight"] == arc_c["lineweight"]
                if case1 and case2  and case3 and (arc_s not in student_lines):
                    correct+=1
                    student_arcs.append(arc_s)
                    arc_flag = True
                    if not case4a or not case4b:
                        correct - ((color_error_penalty/100) * 1)
                        color_error = True
                    if not case5a or not case5b:
                        correct - ((lw_error_penalty/100) * 1)
                        lw_error= True
                    break
        if not arc_flag:
            arc_flags = arc_flags+1
#        if the arc is missing and verbose is 1, print which arc is missing
        if not arc_flag:
            missing_ents = missing_ents +1
        if verbose!=0:
            if not arc_flag:
                print("missing Arc: ",arc_c)
            else:
                if verbose==2:
                    print("found Arc: ",arc_c)
#    Now go through each arc in the student file and see if it was added to the list of correct arcs
#    if there is a arc that was not correct, it will be considered extra and points will be deducted
    for arc_s in entities_student["Arcs"]:
        arc_s = entities_student["Arcs"][arc_s]
        if arc_s not in student_arcs:
            extra_ents = extra_ents +1
            if verbose!=0:
                print("extra Arc: ",arc_s)
    while extra_ents != 0 and arc_flags !=0:
        extra_ents = extra_ents -1
        arc_flags = arc_flags -1

#    see description above for line and arc
    student_circs = []
    for x in range (0, len(entities_correct["Circles"])):
        cir_c = entities_correct["Circles"]["Circle " + str(x)]
        ent_count += 1
        cir_flag = False
        for y in range (0, len(entities_student["Circles"])):
            cir_s = entities_student["Circles"]["Circle " + str(y)]
            if cir_s not in student_circs:
                case1 = round(abs(cir_c["radius"]), thresh) == round(abs(cir_s["radius"]), thresh)
                if len(entities_correct["Lines"])>0 and len(student_lines) >0:
                    pointAc = [cir_c["center"][0], cir_c["center"][1]]
                    pointlc = [entities_correct["Lines"]["Line " + str(0)]["end"][0], entities_correct["Lines"]["Line " + str(0)]["end"][1]]
                    pointAs = [cir_s["center"][0], cir_s["center"][1]]
                    pointls = [student_lines[0]["end"][0], student_lines[0]["end"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case2 = True
                    else:
                        case2 = False
                elif len(entities_correct["Arcs"]) > 0 and len(student_arcs) > 0:
                    z = len(entities_correct["Arcs"])
                    q = len(student_arcs)
                    pointAc = [cir_c["center"][0], cir_c["center"][1]]
                    pointlc = [entities_correct["Arcs"]["Arc " + str(0)]["center"][0], entities_correct["Arcs"]["Arc " + str(0)]["center"][1]]
                    pointAs = [cir_s["center"][0], cir_s["center"][1]]
                    pointls = [student_arcs[0]["center"][0], student_arcs[0]["center"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case2 = True
                    else:
                        case2 = False
                else:
                    case2 = True
            case3a= entities_student["Layers"][cir_s["layer"]]["color"] == entities_correct["Layers"][cir_c["layer"]]["color"]
            case4a = entities_student["Layers"][cir_s["layer"]]["lineweight"] == entities_correct["Layers"][cir_c["layer"]]["lineweight"]
            case3b = cir_s["color"] == cir_c["color"]
            case4b = cir_s["lineweight"] == cir_c["lineweight"]
            if case1 and case2 and (cir_s not in student_circs):
                correct+=1
                student_circs.append(cir_s)
                cir_flag = True
                if not case3a or not case3b:
                    correct - ((color_error_penalty/100) * 1)
                    color_error = True
                if not case4a or not case4b:
                    correct - ((lw_error_penalty/100) * 1)
                    lw_error = True
                break
        if not cir_flag:
            circ_flags = circ_flags+1
        if not cir_flag:
            missing_ents = missing_ents +1
        if verbose!=0:
            if not cir_flag:
                print("missing Circle: ",cir_c)
            else:
                if verbose==2:
                    print("found Circle: ",cir_c)
    for cir_s in entities_student["Circles"]:
        cir_s = entities_student["Circles"][cir_s]
        if cir_s not in student_circs:
            extra_ents = extra_ents +1
            if verbose!=0:
                print("extra Circle: ",cir_s)
    while extra_ents != 0 and circ_flags !=0:
        extra_ents = extra_ents -1
        circ_flags = circ_flags -1

#    see description above for line and arc
    student_ellipses = []
    for x in range(0, len(entities_correct["Ellipses"])):
        ell_c = entities_correct["Ellipses"]["Ellipse " + str(x)]
        ent_count += 1
        ell_flag = False
        for y in range(0, len(entities_student["Ellipses"])):
            ell_s = entities_student["Ellipses"]["Ellipse " + str(y)]
            if ell_s not in student_ellipses:
                case1 = round(math.dist([0,0,0], ell_s["minor axis"])) == round(math.dist([0,0,0], ell_c["minor axis"]))
                case2 = round(math.dist([0,0,0], ell_s["major axis"])) == round(math.dist([0,0,0], ell_c["major axis"]))
                if len(entities_correct["Lines"])>0 and len(student_lines) >0:
                    pointAc = [ell_c["center"][0], ell_c["center"][1]]
                    pointlc = [entities_correct["Lines"]["Line " + str(0)]["end"][0], entities_correct["Lines"]["Line " + str(0)]["end"][1]]
                    pointAs = [ell_s["center"][0], ell_s["center"][1]]
                    pointls = [student_lines[0]["end"][0], student_lines[0]["end"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case3 = True
                    else:
                        case3 = False
                elif len(entities_correct["Arcs"]) > 0 and len(student_arcs) > 0:
                    pointAc = [ell_c["center"][0], ell_c["center"][1]]
                    pointlc = [entities_correct["Arcs"]["Arc " + str(0)]["center"][0], entities_correct["Arcs"]["Arc " + str(0)]["center"][1]]
                    pointAs = [ell_s["center"][0], ell_s["center"][1]]
                    pointls = [student_arcs[0]["center"][0], student_arcs[0]["center"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case3 = True
                    else:
                        case3 = False
                elif len(entities_correct["Circles"]) > 0 and len(student_circs) > 0:
                    pointAc = [cir_c["center"][0], cir_c["center"][1]]
                    pointlc = [entities_correct["Circles"]["Circle " + str(0)]["center"][0], entities_correct["Circles"]["Circle " + str(0)]["center"][1]]
                    pointAs = [cir_s["center"][0], cir_s["center"][1]]
                    pointls = [student_circs[0]["center"][0], student_circs[0]["center"][1]]
                    distc = round(math.dist(pointAc, pointlc),thresh)
                    dists = round(math.dist(pointAs, pointls),thresh)
                    if dists == distc:
                        case3 = True
                    else:
                        case3 = False
                else:
                    case3 = True                     
            case4a = entities_student["Layers"][ell_s["layer"]]["color"] == entities_correct["Layers"][ell_c["layer"]]["color"]
            case5a = entities_student["Layers"][ell_s["layer"]]["lineweight"] == entities_correct["Layers"][ell_c["layer"]]["lineweight"]
            case4b = ell_s["color"] == ell_c["color"]
            case5b = ell_s["lineweight"] == ell_c["lineweight"]
            if case1 and case2 and case3 and (ell_s not in student_ellipses):
                correct+=1
                student_ellipses.append(ell_s)
                ell_flag = True
                if not case4a or not case4b:
                    correct - ((color_error_penalty/100) * 1)
                    color_error = True
                if not case5a or not case5b:
                    correct - ((lw_error_penalty/100) * 1)
                    lw_error = True
                break
        if not ell_flag:
            ell_flags= ell_flags+1
        if not ell_flag:
            missing_ents = missing_ents +1
        if verbose!=0:
            if not ell_flag:
                print("missing Ellipse: ",ell_c)
            else:
                if verbose==2:
                    print("found Ellipse: ",ell_c)
    for ell_s in entities_student["Ellipses"]:
        ell_s = entities_student["Ellipses"][ell_s]
        if ell_s not in student_ellipses:
            extra_ents = extra_ents +1
            if verbose!=0:
                print("extra Ellipse: ",ell_s)
    while extra_ents != 0 and ell_flags !=0:
        extra_ents = extra_ents -1
        ell_flags = ell_flags -1
    

#    see description above for line and arc
    student_hatches = []
    for hatch_c in entities_correct["Hatches"]:
         for path_c in entities_correct["Hatches"][hatch_c]:
            ent_count += 1
            area_c = round(get_area(entities_correct["Hatches"][hatch_c][path_c]), thresh-1)
            type_c = entities_correct["Hatches"][hatch_c][path_c]["type"]
            pattern_c = (entities_correct["Hatches"][hatch_c][path_c]["Pattern"])
            color_c = entities_correct["Hatches"][hatch_c][path_c]["Color"]
            sameColor = True
            samePattern = True
            flag = False
            pat_flag = False
            for hatch_s in entities_student["Hatches"]:
                for path_s in entities_student["Hatches"][hatch_s]:
                    area_s = round((factor**2) * get_area(entities_student["Hatches"][hatch_s][path_s]), thresh-1)
                    type_s = entities_student["Hatches"][hatch_s][path_s]["type"]
                    pattern_s = entities_student["Hatches"][hatch_s][path_s]["Pattern"]
                    pattern_s = str(pattern_s)
                    pattern_c = str(pattern_c)
                    color_s = entities_student["Hatches"][hatch_s][path_s]["Color"]
                    if area_c==area_s and type_c==type_s and color_s == color_c and (path_s not in student_hatches):
                        correct += 1
                        student_hatches.append(path_s)
                        flag = True
                        pat_flag = True
                        if pattern_s != pattern_c:
                            correct - ((hatch_error_penalty/100) * 1)
                            pattern_error = True
                            samePattern = False
                        if color_s != color_c:
                            correct - ((color_error_penalty/100) * 1)
                            color_error = True
                            sameColor = False
                        break
                    
                if flag:
                    break
            if verbose!=0:
                if not samePattern:
                    print("missing pattern: " ,entities_correct["Hatches"][hatch_c][path_c]["Pattern"])
                else:
                    if verbose == 2:
                        print("found pattern: " ,entities_correct["Hatches"][hatch_c][path_c]["Pattern"])
                if not samePattern:
                    print("missing color: " ,entities_correct["Hatches"][hatch_c][path_c]["Color"])
                else:
                    if verbose == 2:
                        print("found color: " ,entities_correct["Hatches"][hatch_c][path_c]["Color"])
                if not pat_flag:
                    print("missing Hatch Area: ",entities_correct["Hatches"][hatch_c][path_c])
                else:
                    if verbose==2:
                        print("found Hatch Area: ",entities_correct["Hatches"][hatch_c][path_c])
                        
    for hatch_s in entities_student["Hatches"]:
        for path_s in entities_student["Hatches"][hatch_s]:
            if path_s not in student_hatches:
                if verbose!=0:
                    print("extra Hatch Area: ",entities_student["Hatches"][hatch_s][path_s])
#    no negative grades, and return a percentage]
    if rotations/ln_count >.50:
        rotation= True
    correct = max(correct,0)
    return [round(((correct/ent_count)*100))- ((extra_ents * extra_ent_penalty)), rotation]

# function to grade two dxf files

def grade(dxf_student, dxf_correct, verbose, thresh, extra_ent_penalty, hatch_error_penalty, color_error_penalty, lw_error_penalty, scaling_error, rot_error):
    possible_angles = []
    maxScore=0
    correctFactor = 0
    score= 0
    trueRot = False
    student_ents = get_ents(dxf_student)
    correct_ents = get_ents(dxf_correct)
    factors = find_factor(student_ents, correct_ents, thresh)
    factors.append(1)
    stu_ents = student_ents
    for factor in factors:
        if factor != 0:
            stu_ents = scale(student_ents, factor)
            score, rot  = (grade_ents(stu_ents, correct_ents, factor, verbose, thresh, extra_ent_penalty, hatch_error_penalty, color_error_penalty, lw_error_penalty))
        if score> maxScore and score<101:
            maxScore = score
            correctFactor = factor
            trueRot = rot
    final_score = maxScore
    if round(correctFactor, thresh-1) !=1:
        final_score = final_score - scaling_error
        if verbose == 2:
            print("Scaling error. File needed to be scaled by a factor of: " + str(correctFactor))
    if trueRot is True:
        final_score = final_score - rot_error
        if verbose == 2:
            print("Rotation error.")

    return max(final_score, 0)
# Read in DXF files
#verify()

bind_ip = "152.3.65.15"
if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
