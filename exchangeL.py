SelectedF=20
import datetime
def datetimenow():
    now = datetime.datetime.now()
    print ("Current date and time : ",now.strftime("%Y-%m-%d %H:%M:%S"))
def parsetext():
    datetimenow()
    print ("Hello Mahendra\n ->> Parse Data Set in List to use this function")
    global List
    global element
    pf=open("trainfinal.txt")
    List = pf.readlines()
    element=len(List)
    #print(List[136488])
    #print(List[136496])
    print ("Number of Instance : ",element)
    #myarray = np.asarray(List)
    #print(List)
    pf.close()

def parsedata():
    datetimenow()
    print ("Hello Mahendra\n ->> Parse Data Set in two dimensional Array to use this function")
    pf = open("trainfinal.txt")
    #wt = open("labelCorrectedDataSet.txt","w")
    global List
    ColRange=SelectedF
    #ColRange=nAttr+1
    List = pf.readlines()
    global element
    element=len(List)
    print (element)
    global a
    '''
    for i in range(len(List)):
        if(i==136488):
            print List[i]
        elif(i==136496):
            print List[i]
        else:
            wt.write(List[i])
    '''
    a=[ [0 for j in range(ColRange)] for i in range(len(List))]
    for i in range(len(List)):
        str=List[i].split(",")
        for j in range(ColRange):
            a[i][j]=str[j]
    List[:]=[]
    print("Length of List : ",len(List))
#    for i in range(element):
#        item=a[i][nAttr]
#        if(i==136488 or i==136496):
#            continue
#        else:
#            wt.write("%s" %item)
    pf.close()
    #wt.close()
    ut = open("rdkddcupU.txt","w")
    for i in range(element):
        for j in range(SelectedF-1):
            ut.write(a[i][j])
            if(j != (SelectedF-2)):
                ut.write(",")
            else:
                ut.write("\n")
    ut.close()
    print("Process completed...")

def checkdatadetail():
    datetimenow()
    print (" Hello Mahendra\n ->> Check data set details to use this function")
    parsetext()
    numberOfClass=3

    p=0
    check=[ 0 for i in range(numberOfClass)]
    totLabelCont=[ 0 for u in range(numberOfClass)]

    nclass=0
    for k in range(element):
        size=p+1
        put=0

        for s in range(size):
            if (check[0]==0):
                check[p]=List[k]

                p=p+1
                put=0
                nclass=nclass+1
                break
            elif check[s]==List[k]:
                put=0
                break
            else:
                str1=List[k]
                put=1
        if put == 1:
            check[p]=str1
            nclass=nclass+1
            p=p+1

        for t in range(numberOfClass):
            if check[t]==List[k]:
                totLabelCont[s]=totLabelCont[s]+1

    for j in range(numberOfClass):
        print(totLabelCont[j],check[j])

    print ("Total Distinct Labels : ",nclass)

def exchangedata():
    datetimenow()
    print (" Hello Mahendra\n ->> Rename labeled attacks to belongs five attacks")
    parsetext()
    for i in range(element):
        x=List[i]
        if (x=="normal.\n"):
            List[i]="normal"

        elif (x=="apache2.\n" or x=="back.\n" or x=="land.\n" or x=="mailbomb.\n" or x=="neptune.\n" or x=="pod.\n" or x=="processtable.\n" or x== "smurf.\n" or x=="teardrop.\n" or x=="udpstorm.\n" ):
             List[i]="DoS"

        elif (x=="ipsweep.\n" or x=="mscan.\n" or x=="nmap.\n" or x=="portsweep.\n" or x=="saint.\n" or x=="satan.\n"):
             List[i]="Probe"

        elif (x=="ftp_write.\n" or x=="guess_passwd.\n" or x=="httptunnel.\n" or x=="imap.\n" or x=="multihop.\n" or x=="named.\n" or x=="phf.\n" or x=="sendmail.\n" or x=="snmpgetattack.\n" or x=="spy.\n" or x=="snmpguess.\n" or x=="warezclient.\n" or x=="warezmaster.\n" or x=="worm.\n" or x=="xlock.\n" or x=="xsnoop.\n"):
             List[i]="R2L"

        elif (x=="buffer_overflow.\n" or x=="loadmodule.\n" or x=="perl.\n" or x=="ps.\n" or x=="rootkit.\n" or x=="sqlattack.\n" or x=="xterm.\n"):
            List[i]="U2R"

        else:
            break
    fp=open("labelCh.txt","w")

    for i in range(element):
        fp.write(List[i])
        fp.write("\n")

    print ("Process Completed...")

def writeintotext():
    parsetext()
    datetimenow()
    print(" Hello Mahendra\n ->> Generate Distinct Data Set")
    fp=open("trainfinal.txt","w")
    #op=open("output.txt","w")
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
            #op.write(element-i)
            print("Remaining for dublicate elimination : ",(element-i))
    fp.write(List[ele])
    fp.close()
    #op.close()
    print("Process Completed...")

#parsetext()
#parsedata()
#writeintotext()
checkdatadetail()
#exchangedata()