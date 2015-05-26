def getDeflections(massmodels,points):
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    if type(massmodels)!=type([]):
        massmodels = [massmodels]
    x0 = x.copy()
    y0 = y.copy()
    for massmodel in massmodels:
        xmap,ymap = massmodel.deflections(x,y)
        y0 -= ymap#.reshape(y.shape)#/scale
        x0 -= xmap#.reshape(x.shape)#/scale
    return x0.reshape(x.shape),y0.reshape(y.shape)



def lens_images(massmodels,sources,points,factor=1,getPix=False,csub=11):
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    if type(massmodels)!=type([]):
        massmodels = [massmodels]
#    x0 = x.flatten()
#    y0 = y.flatten()
    x0 = x.copy()
    y0 = y.copy()
    for massmodel in massmodels:
        xmap,ymap = massmodel.deflections(x,y)
        y0 -= ymap#.reshape(y.shape)#/scale
        x0 -= xmap#.reshape(x.shape)#/scale
    x0,y0 = x0.reshape(x.shape),y0.reshape(y.shape)
   
    if type(sources)!=type([]):
        sources = [sources]
    out = x*0.
    for src in sources:
        out += src.pixeval(x0,y0,factor,csub=csub)

    if getPix==True:
        return out, x0,y0
    return out


def dblPlane(scales,massmodels,sources,points,factor,csub=11):
    if type(points)==type([]):
        x1,y1 = points[0].copy(),points[1].copy()
    else:
        y1,x1 = points[0].copy(),points[1].copy()
    out = x1*0.
    ax_1 = x1*0.
    ay_1 = x1*0.
    for l in massmodels[0]:
        xmap,ymap = l.deflections(x1,y1)
        ax_1 += xmap.reshape(ax_1.shape)
        ay_1 += ymap.reshape(ay_1.shape)
    for s in sources[0]:
        out += s.pixeval(x1-ax_1,y1-ay_1,factor,csub=csub)
    x2 = x1-scales[0,0]*ax_1
    y2 = y1-scales[0,0]*ay_1
    ax_2 = x2*0.
    ay_2 = y2*0.
    for l in massmodels[1]:
        xmap,ymap = l.deflections(x2,y2)
        ax_2 += xmap.reshape(ax_2.shape)
        ay_2 += ymap.reshape(ay_2.shape)
    for s in sources[1]:
        out += s.pixeval(x1-ax_1-ax_2,y1-ay_1-ay_2,factor,csub=csub)
    return out


def multiplePlanes(scales,massmodels,points):
    from numpy import zeros,eye,triu_indices
    if type(points)==type([]):
        x,y = points[0].copy(),points[1].copy()
    else:
        y,x = points[0].copy(),points[1].copy()
    out = x*0.
    nplanes = len(massmodels)
    tmp = scales.copy()
    scales = eye(nplanes)
    scales[triu_indices(nplanes,1)] = tmp
    ax = zeros((x.shape[0],x.shape[1],nplanes))
    ay = ax.copy()
    x0 = x.copy()
    y0 = y.copy()
    out = []
    for p in range(nplanes):
        for massmodel in massmodels[p]:
            xmap,ymap = massmodel.deflections(x0,y0)
            ax[:,:,p] += xmap.reshape(x.shape)
            ay[:,:,p] += ymap.reshape(y.shape)
        x0 = x-(ax[:,:,:p+1]*scales[:p+1,p]).sum(2)
        y0 = y-(ay[:,:,:p+1]*scales[:p+1,p]).sum(2)
        out.append([x0,y0])
    return out
