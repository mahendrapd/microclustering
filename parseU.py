import random
import datetime
import numpy as np
import math
nAttr = 21
SelectedF = 20
setF = math.ceil(SelectedF / 2)
setF = int(setF)
nCluster = 31
nClass = 3
def datetimenow():
    now = datetime.datetime.now()
    print ("Current date and time : ",now.strftime("%Y-%m-%d %H:%M:%S"))

def randnum():
    print ("Hello Mahendra\n ->> Random Number Generation to use this function")
    LR = random.sample(range(0, SelectedF), setF)
    print (LR)
    global L
    L = LR.sort()
    print (LR)
    L = [0 for m in range(setF)]
    j = 0
    k = 0
    for i in range(SelectedF):
        if(k < setF and LR[k] == i):
            k = k + 1
        else:
            L[j] = i
            j = j + 1
    print (L)
    print ("Random process completed..")

def isnumeric(s):
    '''returns True if string s is numeric'''
    return all(c in "0123456789.+-\n" for c in s) and any(c in "0123456789" for c in s)

def parsetext():
    #datetimenow()
    print ("Hello Mahendra\n ->> Parse Data Set in List to use this function")
    global List
    global element
    pf=open("cCindex.txt")
    List = pf.readlines()
    element=len(List)
    #print(List[136488])
    #print(List[136496])
    print ("Number of Instance :",element)
    #myarray = np.asarray(List)
    #print(List[1])
    pf.close()

def countUnClust():
    print ("Parse Data Set in List to use this function")
    global List
    pf=open("subCluster.txt")
    List = pf.readlines()
    element=len(List)
    print ("Number of Instance :",element)
    b = [0 for i in range(element)]
    c = 0
    for i in range(element):
        b[i] = int(List[i])
        if(b[i] == -1):
            c = c + 1
    print("Number of UnClustered element : ", c)
    pf.close()

def parsetextmatchtrain():
    parsetext()
    datetimenow()
    print (" Hello Mahendra\n ->> Parse Original Normalized data and feature selected data for match with train data")
    pm=open("train.txt")
    ptf=open("frequencyOfData.txt","w")
    ListM=pm.readlines()
    elementM=len(ListM)
    #freq=[0 for i in range(element)]
    print ("Number of Instances Original Data Set : ",elementM)
    itn=0
    for i in range(element):
        num=0
        for j in range(elementM):
            if(List[i]==ListM[j]):
                num=num+1
        ptf.write(str(num))
        ptf.write("\n")
        itn=itn+1
        if(itn==2000):
            print ("Iteration Remaining : ",element-i)
            itn=0
    print ("Process Completed..")
    pm.close()
    ptf.close()

def countOutlier():
    #datetimenow()
    print ("Parse Data Set in two dimensional Array and count Outliers")
    pf = open("testfinal.txt")
    global List
    ColRange = SelectedF + 1
    #ColRange=nAttr
    List = pf.readlines()
    global element
    element = len(List)
    print (element)
    global a
    a = [[0 for j in range(ColRange)] for i in range(len(List))]
    for i in range(len(List)):
        strt = List[i].split(",")
        for j in range(ColRange):
            a[i][j] = strt[j]
    List[:] = []
    pf.close()
    #wt.close()
    print("data parsing Process completed...")

def parsedata():
    print ("Hello Mahendra\n ->> Parse Data Set in two dimensional Array to use function")
    pf = open("testfinal.txt")
    global List
    ColRange = SelectedF + 1
    #ColRange=nAttr
    List = pf.readlines()
    global element
    element = len(List)
    print (element)
    global a
    a = [[0 for j in range(ColRange)] for i in range(len(List))]
    for i in range(len(List)):
        strt = List[i].split(",")
        for j in range(ColRange):
            a[i][j] = strt[j]
            if(j < ColRange - 1):
                a[i][j] = float(a[i][j])
    List[:] = []
    pf.close()
    #wt.close()
    print("data parsing Process completed...")

def parseCdata():
    datetimenow()
    print ("Hello Mahendra\n ->> Parse cluster center Data")
    pf = open("clustCenter.txt")
    ColRange=SelectedF
    Listc = pf.readlines()
    numC=len(Listc)
    print (numC)
    global c
    c=[ [0 for j in range(ColRange)] for i in range(numC)]
    for i in range(numC):
        str=Listc[i].split(",")
        for j in range(ColRange):
            c[i][j]=str[j]
    Listc[:]=[]
    pf.close()
    print("Process completed...")

def ConvertS2Udataset():
    datetimenow()
    print ("Hello Mahendra\n ->> Convert Supervised dataset to Unsupervised dataset")
    pf = open("datanorm.txt")
    wt = open("udatanorm.txt","w")
    global List
    labelPosition=41
    ColRange=labelPosition+1
    List = pf.readlines()
    global element
    element=len(List)
    print(element)
    global a
    a=[ [0 for j in range(ColRange)] for i in range(len(List))]
    for i in range(len(List)):
        str=List[i].split(",")
        for j in range(ColRange):
            a[i][j]=str[j]
    List[:]=[]

    for i in range(element):
        for j in range(ColRange-1):
            wt.write(a[i][j])
            if(j != (ColRange-2)):
                wt.write(",")
            else:
                wt.write("\n")
    pf.close()
    wt.close()

def writeintotext():
    datetimenow()
    print ("Hello Mahendra\n ->> Generate Distinct Data Set")
    parsetext()
    fp=open("trainfinal.txt","w")
    ele=element-1
    p=0
    for i in range(ele):
        p=p+1
        notdup=0
        k=i+1
        for j in range(k,element):
            if(List[i]==List[j]):
                notdup=1
                break
        if(notdup == 0):
            fp.write(List[i])
        if(p==5000):
            p=0
            print ("Remaining : ", element-i)
    fp.write(List[ele])
    fp.close()
    print ("Process Completed...")

def traintestclassify():
    datetimenow()
    print ("Hello Mahendra\n ->> Seprate Train and Test Data Sets to use this function")
    #parsedata()
    parsetext()
    wtr=open("train.txt","w")
    wts=open("testfinal.txt","w")
    wti=open("testIndex.txt","w")
    #etest=int(round(float(element)/10))
    etest = 60485
    etrain=element-etest
    dk=random.sample(range(0,element),etest)
    dr=[0 for i in range(etest)]
    for i in range(etest):
        dr[i]=dk[i]
    dk[:]=[]
    '''
    for i in range(etest-1):
        j=i+1
        k=j
        for j in range(etest):
            if(dr[k] > dr[j]):
                dr[k]=dr[j]
                k=j
        if(dr[i] > dr[k]):
            x=dr[i]
            dr[i]=dr[k]
            dr[k]=x
    '''
    dr.sort()

    for i in range(etest):
        wti.write(str(dr[i]))
        wti.write("\n")
        #print dr[i]

    print (etest,etrain,element)

    x=0
    p=0
    for i in range(element):
        if( i != dr[x] ):
            p=p+1
            wtr.write(List[i])
        else:
            wts.write(List[i])
            if(x < etest-1):
                x=x+1
    wtr.close()
    wts.close()

    print ("Process Completed...")

def findelement():
    datetimenow()
    print (" Hello Mahendra\n ->> Check data set attributes details to use this function")
    parsedata()
    labelPosition=21
    point=open("depend.txt","w")
    b=[0 for i in range(element)]
    global ce
    ce=[0 for j in range(labelPosition)]
    t=0
    for j in range(labelPosition):
        p=0
        for i in range(element):
            true=1
            if (p==0):
                b[p]=a[i][j]
                p=p+1
                true=0
            else:
                for k in range(p):
                    if(b[k]==a[i][j]):
                      true=0
                      break
                    else:
                        true=1
            if(true==1):
                b[p]=a[i][j]
                p=p+1

        point.write("sequence id : ")
        point.write(str(j))
        point.write("\n")
        for m in range(p):
            point.write(b[m])
            point.write(" ")
        point.write("\n")
        point.write("Total : ")

        ce[t]=p
        t=t+1

        point.write(str(p))
        point.write("\n")
        print (j,p)

    point.close()

    for j in range(t):
        print j,ce[j]

    print ("Find element process Complete...")

def nordataset():
    parsedata()
    print (" Hello Mahendra\n ->> Normalized data set to use this function")
    wt=open("normdata.txt","w")
    mint=open("minval.txt","w")
    maxt=open("maxval.txt","w")
    for j in range(nAttr-1):
        maxx=0
        miny=90000
        for i in range(element):
            a[i][j] = float(a[i][j])
            if(a[i][j] == -1):
                a[i][j] = 0
            if(maxx < a[i][j]):
                maxx = a[i][j]
            if(miny > a[i][j]):
                miny = a[i][j]
        if(maxx > miny):
            for i in range(element):
                a[i][j] = round((float(a[i][j]-miny)/(maxx-miny)),2)
        mint.write(str(miny))
        mint.write("\n")
        maxt.write(str(maxx))
        maxt.write("\n")
        print("process continue : ", j, miny, maxx)
    for i in range(element):
        for j in range(nAttr):
            wt.write(str(a[i][j]))
            if(j!=(nAttr-1)):
                wt.write(",")
    wt.close()
    print ("nordataset process complete..")

def isConvertString2int(f):
    print(f)
    mSize=200
    exactSize=0
    #fptr.write("sequence number\t")
    #fptr.write(str(f))
    #fptr.write("\n")
    mString= [ 0 for r in range(mSize)]
    mInt=[0 for l in range(mSize)]
    for i in range(element):
        t=0
        if(exactSize==0):
            mString[exactSize]=a[i][f]
            mInt[exactSize]=exactSize
            a[i][f]=exactSize
            exactSize=exactSize+1
            t=1
        else:
            for k in range(exactSize):
                if(a[i][f]==mString[k]):
                    a[i][f]=mInt[k]
                    t=1
        if(t==0):
            mString[exactSize]=a[i][f]
            mInt[exactSize]=exactSize
            a[i][f]=exactSize
            exactSize=exactSize+1
    #for p in range(exactSize):
        #fptr.write(str(mString[p]))
        #if(p != exactSize-1):
            #fptr.write(",")
    #fptr.write("\n")
    #for q in range(exactSize):
        #fptr.write(str(mInt[q]))
        #if(q != exactSize-1):
            #fptr.write(",")
    #fptr.write("\n")

def isParseString():
    parsedata()
    print("Data conversion function string to integer")
    wt=open("data.txt","w")
    for m in range(nAttr-1):
        if(isnumeric(a[0][m]) and m != 8):
            print ("Ok")
        else:
            print ("String in processing of feature : ")
            isConvertString2int(m)
    for i in range(element):
        for j in range(nAttr):
            wt.write(str(a[i][j]))
            if(j != (nAttr - 1)):
                wt.write(",")
    print("Process Completed ...")

def fuzzymembership():
    parsedata()
    print("function of membership calculation")
    wt=open("membshipVal.txt","w")
    nTheta=111
    vMemb=[0 for f in range(nTheta)]
    for j in range(nAttr):
        eSize=0
        nData=[0 for m in range(nTheta)]
        nFreq=[0 for l in range(nTheta)]
        for i in range(element):
            t=0
            if(eSize==0):
                nData[eSize]=a[i][j]
                nFreq[eSize]=nFreq[eSize]+1
                eSize=eSize+1
                t=1
            else:
                for p in range(eSize):
                    if(a[i][j]==nData[p]):
                        nFreq[p]=nFreq[p]+1
                        t=1
            if(t==0):
                nData[eSize]=a[i][j]
                nFreq[eSize]=nFreq[eSize]+1
                eSize=eSize+1
        maxVal=0
        for v in range(eSize):
            if(nFreq[v]>maxVal):
                maxVal=nFreq[v]
        vMemb[j]=float(maxVal)/element
        wt.write(str(vMemb[j]))
        wt.write("\n")
        print ("Processed attribute : ",j)
    print ("Process Completed..")

def featureSelection():
    pt=open("membshipVal.txt")
    wt=open("selectedFeatures.txt","w")
    membData=pt.readlines()
    fSize=len(membData)
    print (fSize)
    p=0
    #print membData
    print ("Feature selection process..")
    for k in range(fSize):
        membData[k]=float(membData[k])
        if(membData[k] < 0.995 and membData[k] > 0.05):
            wt.write(str(k))
            wt.write("\n")
            p=p+1
            print (k)
    print ("Number of Selected Features : ", p)

def finalData():
    parsedata()
    print("Function for stored dataset on selected features")
    fp=open("selectedFeatures.txt")
    wt=open("finaldata.txt","w")
    selF=fp.readlines()
    lenFeature=len(selF)
    for m in range(lenFeature):
        selF[m]=int(selF[m])
    for i in range(element):
        for j in range(lenFeature):
            wt.write(str(a[i][selF[j]]))
            if(j!=(lenFeature-1)):
                wt.write(",")
    print ("Process Completed...")

def setFomation():
    parsedata()
    randnum()
    #L = [2, 3, 4, 5, 6, 9, 12, 16, 17, 19]
    print ("sets density identification for cluster center initialization")
    fp = open("setNum.txt", "w")
    ct = open("setCount.txt", "w")
    arrSet=[0 for k in range(element)]
    setNum = 1
    itr = 0

    for i in range(element):
        if(arrSet[i] == 0):
            setCount = 0
            for k in range(element):
                true = 1
                if(arrSet[k] == 0):
                    for j in range(setF):
                        p = L[j]
                        if(a[i][p]!=a[k][p]):
                            true = 0
                            break
                    if(true == 1):
                        arrSet[k] = setNum
                        setCount = setCount + 1
            ct.write(str(setCount))
            ct.write("\n")
            setNum = setNum + 1
        itr = itr + 1
        fp.write(str(arrSet[i]))
        fp.write("\n")
        if(itr == 1000):
            print ("Remaining process", element-i)
            itr = 0
    print ("Process Completed..")

def kmeanClustering():
    parsedata()
    minClustE = 10
    denom = 2
    fp = open("setNum.txt")
    ct = open("setCount.txt")
    fw = open("clustCenter.txt", "w")
    sc = open("subCluster.txt", "w")
    fi = open("cCindex.txt", "w")
    it = open("niteration.txt", "w")
    sNum = fp.readlines()
    sCount = ct.readlines()
    print("This clustering process of subClusters senquentially")
    setLength = len(sCount)
    for i in range(setLength):
        sCount[i] = int(sCount[i])
    for i in range(element):
        sNum[i] = int(sNum[i])
        for j in range(SelectedF):
            a[i][j]=float(a[i][j])
    print(len(sNum), len(sCount))
    fSubClust = [-1 for i in range(element)]
    #finalCluster=[-1 for i in range(element)]
    clusteredElement = 0
    sLength = -1
    loopTrue = 1
    while(loopTrue == 1):
        sLength = sLength + 1
        if(sLength >= setLength):
            break
        if(sCount[sLength] >= minClustE):
            avgDistSC = 0
            cCenter = [0 for c in range(SelectedF)]
            ecludDist = [0 for i in range(element)]
            includedE = [-1 for i in range(element)]
            for i in range(element):
                if(sNum[i] == sNum[sLength]):
                    #includedE[i]=sNum[sLength]
                    for j in range(SelectedF):
                        cCenter[j] = cCenter[j] + a[i][j]
            for j in range(SelectedF):
                cCenter[j] = float(cCenter[j]) / sCount[sLength]
            navg = 0
            for i in range(element):
                dist = 0
                for j in range(SelectedF):
                    dist = dist + (a[i][j] - cCenter[j])**2
                ecludDist[i] = math.sqrt(dist)
                if(sNum[i] == sNum[sLength]):
                    navg = navg + 1
                    avgDistSC = avgDistSC + ecludDist[i]
            avgDistSC = float(avgDistSC) / (denom * navg)
            for i in range(element):
                if(ecludDist[i] <= avgDistSC):
                    includedE[i] = sNum[sLength]
            true = 1
            iteration = 0
            while(true == 1):
                cCenter = [0 for c in range(SelectedF)]
                incClustE=[-1 for i in range(element)]
                iteration = iteration + 1
                subCount = 0
                nSubClust = 0
                for i in range(element):
                    if(includedE[i] == sNum[sLength]):
                        subCount = subCount + 1
                        for j in range(SelectedF):
                            cCenter[j] = cCenter[j] + a[i][j]
                if(subCount == 0):
                    break
                for j in range(SelectedF):
                    cCenter[j] = float(cCenter[j]) / subCount
                for i in range(element):
                    dist = 0
                    for j in range(SelectedF):
                        dist = dist+(a[i][j] - cCenter[j])**2
                    ecludDist[i] = math.sqrt(dist)
                    if(ecludDist[i] <= avgDistSC and fSubClust[i] == -1):
                        nSubClust = nSubClust + 1
                        incClustE[i] = sNum[sLength]
                loopTest = 0
                for i in range(element):
                    if(incClustE[i] != includedE[i]):
                       loopTest = loopTest + 1
                       includedE[i] = incClustE[i]
                if(loopTest == 0):
                    clusteredElement = clusteredElement + nSubClust
                    it.write(str(iteration))
                    it.write("\n")
                    true = 0
                    for i in range(element):
                        if(fSubClust[i] == -1 and includedE[i] != -1):
                            fSubClust[i] = includedE[i]
            if(nSubClust > 2):
                fi.write(str(sNum[sLength]))
                fi.write("\n")
                for j in range(SelectedF):
                    cCenter[j] = round(cCenter[j], 2)
                    fw.write(str(cCenter[j]))
                    if(j != (SelectedF - 1)):
                        fw.write(",")
                fw.write("\n")
                print("Procecss : ", (element-clusteredElement), nSubClust, clusteredElement, sLength, sNum[sLength])
        #print("Remaining Procecss : ", (element-clusteredElement), nSubClust, clusteredElement, sLength)
        if(clusteredElement >= element or sLength >= setLength):
            loopTrue = 0
    for i in range(element):
        sc.write(str(fSubClust[i]))
        sc.write("\n")
    fp.close()
    ct.close()
    fw.close()
    fi.close()
    sc.close()
    print("Process Completed..")

def minSpanningTree():
    #parsedata()
    fc=open("clustCenter.txt")
    fi=open("cCindex.txt")
    rci=open("cCindexFReduce.txt","w")
    sc=open("subCluster.txt")
    fp=open("fassClust.txt","w")
    finalClust=sc.readlines()
    fcCi=fi.readlines()
    #fcCi=np.asarray(fcCi)
    cList=fc.readlines()
    lenCC=len(fcCi)
    nData=len(finalClust)
    fClust=[0 for i in range(nData)]
    cCi=[0 for k in range(lenCC)]
    c=[ [0 for j in range(SelectedF)] for i in range(lenCC)]
    for i in range(lenCC):
        str=cList[i].split(",")
        for j in range(SelectedF):
            c[i][j]=str[j]
    for i in range(lenCC):
        fcCi[i]=int(fcCi[i])
        cCi[i]=fcCi[i]
        for j in range(SelectedF):
            c[i][j]=float(c[i][j])
    for i in range(nData):
        finalClust[i]=int(finalClust[i])
        fClust[i]=finalClust[i]
    distCC=[[0 for j in range(lenCC)] for i in range(lenCC)]
    minDistCC=[0 for i in range(lenCC)]
    minIndex=[-1 for i in range(lenCC)]
    for k in range(lenCC):
        for i in range(lenCC):
            dist=0
            for j in range(SelectedF):
                dist=dist+(c[k][j]-c[i][j])**2
            distCC[k][i]=math.sqrt(dist)
#    for i in range(lenCC):
#        print("Sequence number : ",i)
#        for j in range(lenCC):
#            print(j,distCC[i][j])
    print("\n MinDistanceIndex for assigning\n")
    for i in range(lenCC):
        minD=10000
        index=-1
        for k in range(lenCC):
            if(i!=k and minD>distCC[i][k]):
                minD=distCC[i][k]
                index=k
        minIndex[i]=index
        minDistCC[i]=distCC[i][index]
        print(i,index,minDistCC[i])
    conD=lenCC
    vis=[0 for i in range(lenCC)]
    while(conD>=1):
        minD=10000
        index=-1
        for i in range(lenCC):
            if(minD>minDistCC[i] and vis[i]==0):
                minD=minDistCC[i]
                index=i
        vis[index]=1
        #conD=conD-1
        true=0
        for k in range(nData):
            if(fClust[k]==cCi[index]):
                true=1
                fClust[k]=cCi[minIndex[index]]
        if(true==1):
            conD=conD-1
            cCi[index]=cCi[minIndex[index]]
                #fcCi[index]=fcCi[minIndex[index]]
    Range=lenCC+1
    cClustCount=[0 for k in range(Range)]
    cClust=[0 for k in range(Range)]
    uBound=0
    for i in range(nData):
        tr=0
        if(uBound==0):
            cClust[uBound]=fClust[i]
            cClustCount[uBound]=cClustCount[uBound]+1
            uBound=uBound+1
            tr=1
        else:
            for k in range(uBound):
                if(fClust[i]==cClust[k]):
                    cClustCount[k]=cClustCount[k]+1
                    tr=1
                    break
        if(tr==0):
            cClust[uBound]=fClust[i]
            cClustCount[uBound]=cClustCount[uBound]+1
            uBound=uBound+1
    print("\n Final Data\n")
    for k in range(Range):
        print(k,cClustCount[k],cClust[k])
    print("\nClusters are ..\n")
    print(fcCi)
    print(cCi)
#................ list write into text file ...........................
    for item in cCi:
        rci.write("%s\n" % item)
    for item in fClust:
        fp.write("%s\n" % item)
#.....................................................................
    fc.close()
    fi.close()
    sc.close()
    fp.close()
    rci.close()
    print("Process Completed..")

def workArrange():
    parsedata()
    parseCdata()
    ci=open("cCindex.txt")
    #rci=open("cCindexFReduce.txt")
    #ti=open("testIndex.txt")
    fp=open("testData2Class.txt","w")
    print("Function for work arrange")
    cCi=ci.readlines()
    lenCi=len(cCi)
    print(lenCi)
    #rcCi=rci.readlines()
    #tIndex=ti.readlines()
    for i in range(element):
        #tIndex[i]=int(tIndex[i])
        for j in range(SelectedF):
            a[i][j]=float(a[i][j])
    for k in range(lenCi):
        cCi[k]=int(cCi[k])
        #rcCi[k]=int(rcCi[k])
        for j in range(SelectedF):
            c[k][j]=float(c[k][j])
    ecludDt=[0 for i in range(lenCi)]
    for k in range(element):
        for i in range(lenCi):
            dist=0
            for j in range(SelectedF):
                dist=dist+(a[k][j]-c[i][j])**2
            ecludDt[i]=math.sqrt(dist)
        tdata=100000
        index=-1
        for i in range(lenCi):
            if(tdata>ecludDt[i]):
                tdata=ecludDt[i]
                index=i
        fp.write(str(cCi[index]))
        fp.write("\n")
    print("Process completed..")
    ci.close()
    #rci.close()
    #ti.close()
    fp.close()

def workTest():
    parsedata()
    print("Arrange work in Confusion Matrix")
    ci = open("cCindex.txt")
    tc = open("testData2Class.txt")
    fr = open("finalResult2.txt", "w")
    cCind = ci.readlines()
    testedClass = tc.readlines()
    #dataLabel=lch.readlines()
    #testDataInd=ti.readlines()
    nData = len(testedClass)
    nClust = len(cCind)
    for k in range(nClust):
        cCind[k] = int(cCind[k])
    for i in range(nData):
        #testDataInd[i] = int(testDataInd[i])
        testedClass[i] = int(testedClass[i])
    print("Test process started ...")
    #print(len(dataLabel))
    testme = 0
    print(nData, nClust)
    size = nCluster
    while(testme < nClust):
        arrLabel = [0 for p in range(size)]
        totCount = [0 for q in range(size)]
        fr.write("\n Sub-cluster : ")
        fr.write(str(cCind[testme]))
        fr.write("\n")
        rB = 0
        for i in range(nData):
            true = 1
            index = -1
            enter = 0
            if(testedClass[i] == cCind[testme]):
                enter = 1
                if(rB == 0):
                    arrLabel[rB] = a[i][SelectedF]
                    totCount[rB] = totCount[rB] + 1
                    rB = rB + 1
                    true = 0
                else:
                    for k in range(rB):
                        if(arrLabel[k] == a[i][SelectedF]):
                            index = k
                            true = 0
                            break
            if(true == 1 and enter == 1):
                arrLabel[rB] = a[i][SelectedF]
                totCount[rB] = totCount[rB] + 1
                rB = rB + 1
            if(index != -1 and enter == 1):
                totCount[index] = totCount[index] + 1
        for m in range(size):
            if(totCount[m] != 0):
                fr.write(str(totCount[m]))
                fr.write("\t")
                fr.write(str(arrLabel[m]))
        testme = testme + 1
    ci.close()
    tc.close()
    #lch.close()
    #ti.close()
    fr.close()
    print("Process Completed...")

def parsecenter():
    print ("Parse Cluster Center Data into two dimensional Array")
    cf = open("clustCenter.txt")
    global Listc
    Listc = cf.readlines()
    print(len(Listc))
    global c
    c = [[0 for j in range(SelectedF)] for i in range(nCluster)]
    for i in range(nCluster):
        strt = Listc[i].split(",")
        for j in range(SelectedF):
            c[i][j] = strt[j]
            c[i][j] = float(c[i][j])
    Listc[:] = []
    cf.close()
    print("data parsing Process completed...")

def testkmeans():
    parsedata()
    parsecenter()
    rt = open("results.txt", "w")
    lt = open("labels.txt", "w")
    result = [[0 for j in range(nClass)] for i in range(nCluster)]
    label = [0 for r in range(nClass)]
    l = 0
    for i in range(element):
        index = 0
        ld = 0
        t = 0
        dist = 99999
        ecluD = [0 for r in range(nCluster)]
        for k in range(nCluster):
            for j in range(SelectedF):
                ecluD[k] = ecluD[k] + (a[i][j] - c[k][j])**2
            ecluD[k] = math.sqrt(ecluD[k])
            if(dist > ecluD[k]):
                dist = ecluD[k]
                index = k
        if(l == 0):
            label[l] = a[i][SelectedF]
            #result[index][l] += 1
            ld = l
            l = l + 1
            t = 1
        else:
            for s in range(nClass):
                if(a[i][SelectedF] == label[s]):
                    ld = s
                    t = 1
                    break
        if(t == 0):
            label[l] = a[i][SelectedF]
            ld = l
            l = l + 1
        result[index][ld] += 1
    print(label)
    for i in range(nClass):
        lt.write(str(label[i]))
    print(result)
    for i in range(nCluster):
        for j in range(nClass):
            rt.write(str(result[i][j]))
            if(j < nClass - 1):
                rt.write("\t")
            else:
                rt.write("\n")
    print("Process Completed ..")

#randnum()
#parsetext()
#parsetextmatchtrain()
#parsedata()
#ConvertS2Udataset()
#findelement()
#nordataset()
#isParseString()
#fuzzymembership()
#featureSelection()
#finalData()
#traintestclassify()
#writeintotext()
#parsetextmatchtrain()
#setFomation()
#kmeanClustering()
#minSpanningTree()
#workArrange()
#workTest()
#testkmeans()
#countUnClust()