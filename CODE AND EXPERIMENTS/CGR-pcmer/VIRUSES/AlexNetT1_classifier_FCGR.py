## MODULE
from module import *

## FUNCTIONS
from functions_Net import create_AlexNet, preprocessing, metrics, saveModel, plot_loss_accuracy, saveConfMatrixClassReport

if __name__ == '__main__':
      
      ## setting parameters
      #type_encoder = "Grayscale"
      type_encoder = "RGB"

      dataset = 'Viruses/Coronaviruses'
      #dataset = 'Viruses/HIV1'
      #dataset = 'Viruses/HIV2'
      #dataset = 'Viruses/Dengue'
      #dataset = 'Viruses/HepatitisC'
      #dataset = 'Viruses/HepatitisB1'
      #dataset = 'Viruses/HepatitisB2'
      #dataset = 'Viruses/InfluenzaA'
     
      k = 4
      threshold = 1
      type_encodingColour = "all encoding MAXFREQUENCYFCGR"
      job = 'FCGR RGB-PCMER-FREQ OPTIMIZATION (k=' + str(k) + ' threshold=' + str(threshold) + " " + type_encodingColour + ')'

      
      X_data, y_data, nb_classes=preprocessing('AlexNet', type_encoder, dataset, k, job)
      batch_size=30
      epoch=30 # 60 for Influenza A
      X_data= np.array(X_data)
      y_data = np.array(y_data)
      print(y_data.shape)
      print(X_data.shape)     
      X_data = X_data.reshape((-1, 227, 227, 3))      
      X_data = X_data.astype('float32')
      print('data shape: {}'.format(X_data.shape))
      print('data labels shape: {}'.format(y_data.shape))
      print('nb_classes: {}'.format(nb_classes))
      shape = X_data.shape[1:]
      print('shape ', shape)

      X = X_data; y = y_data
      # Indexing for training and validation sets (70 + 20)
      # Indexing for testing set (10)
      X_data, X_test, y_data, y_test = train_test_split(X, y, test_size=0.10, random_state=42, stratify = y)  

      skf_AlexNet = StratifiedKFold(n_splits=5,shuffle=True,random_state=20)
      tmp=1; model_AlexNet = Sequential()      

      start = time.time()
      n_task = 2
      for train_index, validation_index in skf_AlexNet.split(X_data, y_data):
          model_AlexNet = Sequential()
          X_train, X_validation = X_data[train_index], X_data[validation_index]
          y_train, y_validation = y_data[train_index], y_data[validation_index]
          print('Fold'+str(tmp)+':')   

          with ProcessPoolExecutor(n_task) as e:
            e.map(create_AlexNet, range(n_task))      
            model_AlexNet = create_AlexNet(model_AlexNet, (227,227,3), nb_classes)
            history=model_AlexNet.fit(X_train[:], y_train[:],
                  batch_size=batch_size,
                  epochs=epoch,
                  validation_data=(X_validation[:], y_validation[:]))
            print('Fold'+str(tmp)+'is finished')
      end = time.time() 

      val_acc = "Validation accuracy" + str((history.history['val_accuracy'])[-1])

      training_time  = "model training time of AlexNet Model with " + type_encoder + " encoder unit: " + str(end-start) + ' s'
      print(training_time)    
          
       # save the classification model as a file
      saveModel(dataset, model_AlexNet, 'AlexNet', type_encoder, job, batch_size, epoch)    

      acc_test = plot_loss_accuracy(history, model_AlexNet, X_test, y_test, dataset, 'AlexNet', job, batch_size, epoch)

      conf_matrix, class_report = metrics(X_test, y_test, model_AlexNet)
      print('\n', conf_matrix, '\n', class_report)

      # save the results of classification model
      saveConfMatrixClassReport('AlexNet', training_time, acc_test, conf_matrix, class_report, dataset, type_encoder, job, val_acc, batch_size, epoch)
