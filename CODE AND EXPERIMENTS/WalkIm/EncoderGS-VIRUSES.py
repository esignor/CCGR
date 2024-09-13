from module import *

def Encode_Grayscale(masterpath, directory):
    mypath=masterpath

    w, h = 1024,1024
    w2, h2 = 64,64
    my_list = os.listdir(mypath)

    for i in range(len(my_list)):
        Subpath=mypath+'/'+my_list[i]
        print(Subpath)
        imgpath=my_list[i]
        onlyfiles = [f for f in listdir(Subpath) if isfile(join(Subpath, f))]
        for s in range(0,len(onlyfiles)):

            data = np.zeros((h, w), dtype=np.uint8)

            x=''
            with open(Subpath + '/' + onlyfiles[s]) as f:
                content = f.readlines()

            content = [x.strip() for x in content]
            for k in range(1,len(content)):
                x+=content[k]
            i=int(w/2)-1 # x
            j=int(h/2)-1 # y

            for k in range (0,len(x)):
                if x[k]=='A' or x[k]=='a': # A (-1,1)
                    i-=1
                    j+=1
                    if(i==-1 or i==w or j==-1 or j==h):
                        i=int((w/2)-1)
                        j=int((h/2)-1)
                    data[i,j]+=255
                if x[k]=='C' or x[k]=='c': # C (-1,-1)
                    i-=1
                    j-=1
                    if(i==-1 or i==w or j==-1 or j==h):
                        i=int((w/2)-1)
                        j=int((h/2)-1)
                    data[i,j]+=255
                if x[k]=='G' or x[k]=='g': # G (1,-1)
                    i+=1
                    j-=1
                    if(i==-1 or i==w or j==-1 or j==h):
                        i=int((w/2)-1)
                        j=int((h/2)-1)
                    data[i,j]+=255
                if x[k]=='T' or x[k]=='t': # T (1,1)
                    i+=1
                    j+=1
                    if(i==-1 or i==w or j==-1 or j==h):
                        i=int((w/2)-1)
                        j=int((h/2)-1)
                        data[i,j]+=255


            img = Image.fromarray(data)
            new_image = img.resize((w2,h2))


            if not os.path.exists(directory+'/'+imgpath):
                os.makedirs(directory+'/'+imgpath)

            imagepath=directory+'/'+imgpath+'/'+onlyfiles[s]+'.png'
            new_image.save(imagepath)



dataset = 'Coronaviruses'
#dataset = 'HIV1'
#dataset = 'HIV2'
#dataset = 'Dengue'
#dataset = 'HepatitisC'
#dataset = 'HepatitisB1'
#dataset = 'HepatitisB2'
#dataset = 'InfluenzaA'

MAIN_FOLDER = 'CODE AND EXPERIMENTS/'

if dataset == 'Coronaviruses':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_Coronaviruses_GrayscaleImg'; masterpath=MAIN_FOLDER + 'DATASET/7classes_coronaviruses dataset'

elif dataset == 'HIV1':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_HIV1_GrayscaleImg'; masterpath=MAIN_FOLDER + 'DATASET/12classes_hiv dataset'

elif dataset == 'HIV2':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_HIV2_GrayscaleImg'; masterpath=MAIN_FOLDER + 'DATASET/37classes_hiv dataset'

elif dataset == 'Dengue':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_Dengue_GrayscaleImg';  masterpath=MAIN_FOLDER + 'DATASET/4classes_dengue dataset'

elif dataset == 'HepatitisC':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_HepatitisC_GrayscaleImg';  masterpath=MAIN_FOLDER + 'DATASET/6classes_hepatitisC dataset'

elif dataset == 'HepatitisB1':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_HepatitisB1_GrayscaleImg';  masterpath=MAIN_FOLDER + 'DATASET/8classes_hepatitisB dataset'

elif dataset == 'HepatitisB2':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_HepatitisB2_GrayscaleImg';  masterpath=MAIN_FOLDER + 'DATASET/13classes_hepatitisB dataset'
    
elif dataset == 'InfluenzaA':
    directory=MAIN_FOLDER + 'WalkIm/OUTWalkIm_InfluenzaA_GrayscaleImg';  masterpath=MAIN_FOLDER + 'DATASET/56classes_influenzaA_reduced dataset'

if not os.path.exists(directory):
   os.makedirs(directory)

if __name__ == '__main__':
    Encode_Grayscale(masterpath, directory)
