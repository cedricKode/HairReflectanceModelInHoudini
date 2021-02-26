# Read Converted Model

node = hou.pwd()
geo = node.geometry()

fileToRead = open(hou.evalParm("hairfile"), "r")
for x in fileToRead:
    strandSplitted = x.split()
    index = 0
    poly = geo.createPolygon()
    while index < len(strandSplitted) - 2:
        pt0 = geo.createPoint()
        pt0.setPosition(hou.Vector3(float(strandSplitted[index]), 
                                    float(strandSplitted[index+1]), 
                                    float(strandSplitted[index+2])))
        poly.addVertex(pt0)
        index += 3
    poly.setIsClosed(False)
fileToRead.close()  