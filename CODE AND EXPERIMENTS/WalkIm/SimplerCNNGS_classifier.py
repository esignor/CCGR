## MODULE
from module import *
## FUNCTIONS
from functions import create_SimplerCNN, preprocessing, metrics, plot_loss_accuracy, saveConfMatrixClassReport, saveModel

if __name__ == '__main__':

  num = 123
  np.random.seed(num)
  os.environ['PYTHONHASHSEED'] = str(num)
  os.environ['TF_DETERMINISTIC_OPS'] = '1'
  os.environ['TF_CUDNN_DETERMINISTIC'] = '1'
  tf.random.set_seed(num)
  tf.keras.utils.set_random_seed(num)
  tf.config.experimental.enable_op_determinism()


  ## setting parameters
  type_encoder = "Grayscale"
  #type_encoder = "RGB"
  dataset = 'Coronaviruses'
  #dataset = 'HIV1'
  #dataset = 'HIV2'
  #dataset = 'Dengue'
  #dataset = 'HepatitisC'
  #dataset = 'HepatitisB1'
  #dataset = 'HepatitisB2'
  #dataset = 'InfluenzaA'
  
  X_data, y_data, nb_classes=preprocessing('CNN', type_encoder, dataset)
  batch_size=64
  epoch=10 # 30 for Influenza A
  X_data= np.array(X_data)
  y_data = np.array(y_data)

  if type_encoder == "Grayscale":
    X_data = X_data.reshape((-1, 64, 32, 2)) # redefine the shape based on the DATASET: (n.examples, _, _, _) must have 4 dimension
  else:
    X_data = X_data.reshape((-1, 64, 64, 3))

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
  
  skf_cnn = StratifiedKFold(n_splits=5,shuffle=True,random_state=20)
  tmp=1; model_cnn = Sequential()

  start = time.time()
  for train_index, validation_index in skf_cnn.split(X_data, y_data):
      model_cnn = Sequential()
      X_train, X_validation = X_data[train_index], X_data[validation_index]
      y_train, y_validation = y_data[train_index], y_data[validation_index]
      print('Fold'+str(tmp)+':')

      model_cnn = create_SimplerCNN(model_cnn, shape, nb_classes)
      model_cnn.summary()
      history=model_cnn.fit(X_train[:], y_train[:],
            batch_size=batch_size,
            epochs=epoch,
            validation_data=(X_validation[:], y_validation[:]))
      print('Fold'+str(tmp)+'is finished')
  end = time.time()

  training_time  = "model training time of Simpler CNN Model with " + type_encoder + " encoder unit: " + str(end-start) + ' s'
  print(training_time)

  # save the classification model as a file
  saveModel(dataset, model_cnn, 'Simpler CNN', type_encoder)


  acc = plot_loss_accuracy(history, model_cnn, X_test, y_test, dataset, 'Simpler CNN', type_encoder)

  conf_matrix, class_report = metrics(X_test, y_test, model_cnn)
  print('\n', conf_matrix, '\n', class_report)

  # save the results of classification model
  saveConfMatrixClassReport('Simpler CNN', training_time, acc, conf_matrix, class_report, dataset, type_encoder)
