import sys
sys.path.insert(1, 'CODE AND EXPERIMENTS/CGR-pcmer/')
import VIRUSES
from VIRUSES.module import *

def create_AlexNet(model, shape, nb_classes):
  #model = Sequential()
  model.add(layers.Conv2D(filters=96, kernel_size=(11, 11), strides=(4, 4), activation="relu", input_shape= shape))
  model.add(layers.BatchNormalization())
  model.add(layers.MaxPool2D(pool_size=(3, 3), strides= (2, 2)))
  model.add(layers.Conv2D(filters=256, kernel_size=(5, 5), strides=(1, 1), activation="relu", padding="same"))
  model.add(layers.BatchNormalization())
  model.add(layers.MaxPool2D(pool_size=(3, 3), strides=(2, 2)))
  model.add(layers.Conv2D(filters=384, kernel_size=(3, 3), strides=(1, 1), activation="relu", padding="same"))
  model.add(layers.BatchNormalization())
  model.add(layers.Conv2D(filters=384, kernel_size=(3, 3), strides=(1, 1), activation="relu", padding="same"))
  model.add(layers.BatchNormalization())
  model.add(layers.Conv2D(filters=256, kernel_size=(3, 3), strides=(1, 1), activation="relu", padding="same"))
  model.add(layers.BatchNormalization())
  model.add(layers.MaxPool2D(pool_size=(3, 3), strides=(2, 2)))
  model.add(layers.Flatten())
  model.add(layers.Dense(4096, activation="relu"))
  model.add(layers.Dropout(0.5, seed=0))
  model.add(layers.Dense(4096, activation="relu"))
  model.add(layers.Dropout(0.5, seed=0))
  model.add(layers.Dense(nb_classes, kernel_regularizer=keras.regularizers.l1(0.01), activation="softmax"))
  model.compile(loss='sparse_categorical_crossentropy', optimizer=tf.optimizers.SGD(learning_rate=0.001), metrics=['accuracy'])
  return model


def convolutional_block(X, f, filters, stage, block, s = 2):
  conv_name_base = 'res' + str(stage) + block + '_branch'
  bn_name_base = 'bn' + str(stage) + block + '_branch'
  F1, F2, F3 = filters
  X_shortcut = X
  X = layers.Conv2D(filters = F1, kernel_size = (1, 1), strides = (s,s), name = conv_name_base + '2a', kernel_initializer = glorot_uniform(seed=0))(X)
  X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2a')(X)
  X = layers.Activation('relu')(X)
  X = layers.Conv2D(filters = F2, kernel_size = (f, f), strides = (1,1), padding = 'same', name = conv_name_base + '2b', kernel_initializer = glorot_uniform(seed=0))(X)
  X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2b')(X)
  X = layers.Activation('relu')(X)
  X = layers.Conv2D(filters = F3, kernel_size = (1, 1), strides = (1,1), padding = 'valid', name = conv_name_base + '2c', kernel_initializer = glorot_uniform(seed=0))(X)
  X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2c')(X)
  X_shortcut = layers.Conv2D(filters = F3, kernel_size = (1, 1), strides = (s,s), padding = 'valid', name = conv_name_base + '1', kernel_initializer = glorot_uniform(seed=0))(X_shortcut)
  X_shortcut = layers.BatchNormalization(axis = 3, name = bn_name_base + '1')(X_shortcut)
  X = layers.Add()([X_shortcut, X])
  X = layers.Activation('relu')(X)

  return X


def identity_block(X, f, filters, stage, block):

    conv_name_base = 'res' + str(stage) + block + '_branch'
    bn_name_base = 'bn' + str(stage) + block + '_branch'

    F1, F2, F3 = filters

    X_shortcut = X

    X = layers.Conv2D(filters = F1, kernel_size = (1, 1), strides = (1,1), padding = 'valid', name = conv_name_base + '2a', kernel_initializer = glorot_uniform(seed=0))(X)
    X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2a')(X)
    X = layers.Activation('relu')(X)

    X = layers.Conv2D(filters = F2, kernel_size = (f, f), strides = (1,1), padding = 'same', name = conv_name_base + '2b', kernel_initializer = glorot_uniform(seed=0))(X)
    X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2b')(X)
    X = layers.Activation('relu')(X)

    X = layers.Conv2D(filters = F3, kernel_size = (1, 1), strides = (1,1), padding = 'valid', name = conv_name_base + '2c', kernel_initializer = glorot_uniform(seed=0))(X)
    X = layers.BatchNormalization(axis = 3, name = bn_name_base + '2c')(X)

    # Add shortcut value to main path
    X = layers.Add()([X_shortcut, X])
    X = layers.Activation('relu')(X)

    return X

def create_ResNet50(model, shape, nb_classes):
    X_input = layers.Input(shape)
    X = layers.ZeroPadding2D((3, 3))(X_input)
    X = layers.Conv2D(64, (7, 7), strides = (2, 2), name = 'conv1', kernel_initializer = glorot_uniform(seed=0))(X)
    X = layers.BatchNormalization(axis = 3, name = 'bn_conv1')(X)
    X = layers.Activation('relu')(X)
    X = layers.MaxPooling2D((3, 3), strides=(2, 2))(X)
    X = convolutional_block(X, f = 3, filters = [64, 64, 256], stage = 2, block='a', s = 1)
    X = identity_block(X, 3, [64, 64, 256], stage=2, block='b')
    X = identity_block(X, 3, [64, 64, 256], stage=2, block='c')
    X = convolutional_block(X, f = 3, filters = [128, 128, 512], stage = 3, block='a', s = 2)
    X = identity_block(X, 3, [128, 128, 512], stage=3, block='b')
    X = identity_block(X, 3, [128, 128, 512], stage=3, block='c')
    X = identity_block(X, 3, [128, 128, 512], stage=3, block='d')
    X = convolutional_block(X, f = 3, filters = [256, 256, 1024], stage = 4, block='a', s = 2)
    X = identity_block(X, 3, [256, 256, 1024], stage=4, block='b')
    X = identity_block(X, 3, [256, 256, 1024], stage=4, block='c')
    X = identity_block(X, 3, [256, 256, 1024], stage=4, block='d')
    X = identity_block(X, 3, [256, 256, 1024], stage=4, block='e')
    X = identity_block(X, 3, [256, 256, 1024], stage=4, block='f')
    X = convolutional_block(X, f = 3, filters = [512, 512, 2048], stage = 5, block='a', s = 2)
    X = identity_block(X, 3, [512, 512, 2048], stage=5, block='b')
    X = identity_block(X, 3, [512, 512, 2048], stage=5, block='c')
    X = layers.AveragePooling2D(pool_size=(2, 2),name='avg_pool')(X)
    X = layers.Flatten()(X)
    X = layers.Dense(nb_classes, kernel_regularizer=keras.regularizers.l1(0.01), activation='softmax', name='fc' + str(nb_classes), kernel_initializer = glorot_uniform(seed=0))(X)
    model = Model(inputs = X_input, outputs = X, name='ResNet50')
    model.compile(loss='sparse_categorical_crossentropy', optimizer=tf.optimizers.SGD(learning_rate=0.001), metrics=['accuracy'])
    return model

  
def preprocessing(type_arch, type_encoder, dataset, k, job):
    MAIN_FOLDER = 'CODE AND EXPERIMENTS/'
    out_directory = MAIN_FOLDER + 'CGR-pcmer/VIRUSES/'

## for the resize method bicubic resampling is used (default value): there is a fixed
# number of source image pixels for each target pixel (i.e., 4*4 for bicubic) where 
# affine transformations are applied.
    X = []
    y = []
    path_main = out_directory + 'OUTFCGRPCMER_ENCODER/' + dataset + '/' + job
    os.chdir(path_main)
    dirs = filter(os.path.isdir, os.listdir(os.curdir))
    for dir in dirs:
        path_to_subdir = str(dir)
        for im_path in os.listdir(path_to_subdir):
            im_frame = Image.open(path_to_subdir + '/' + im_path)

            if type_arch == "AlexNet":
              ## resize image for AlexNet
              single_channel = im_frame.resize((227, 227))
              # np_frame = np.array((im_frame.resize((227, 227)).getdata()))
              if type_encoder == "RGB":
                np_frame = np.concatenate([single_channel], axis = -1)
              else:
                np_frame = np.concatenate([single_channel, single_channel, single_channel], axis = -1)


            elif type_arch == "ResNet50":
              ## resize image for ResNet50
              single_channel = im_frame.resize((224, 224))
              # np_frame = np.array((im_frame.resize((224, 224)).getdata()))
              if type_encoder == "RGB":
                np_frame = np.concatenate([single_channel], axis = -1)
              else:
                np_frame = np.concatenate([single_channel, single_channel, single_channel], axis = -1)

            X.append(np_frame)
            y.append(path_to_subdir.split('.')[0])
    
    os.chdir('../../../../')
    unique = list(dict.fromkeys(y))
    dct = {}
    cnt = 0
    for lab in unique:
        dct[str(lab)] = cnt
        cnt += 1

    nb_classes = len(dct)
    new_label = []
    for l in y:
        new_label.append(dct[l])

    y = new_label
    return X, y, nb_classes
    
    
    
def plot_loss(history, dataset, model_net, job, batch_size, epoch):
  directory= 'CGR-pcmer ' + dataset.split('/')[1] + ' Results/'
  if not os.path.exists(directory):
    os.makedirs(directory)

  plt.figure(figsize=(10,6))
  plt.plot(history.epoch,history.history['loss'], label = "Training loss")
  plt.plot(history.epoch,history.history['val_loss'], label = "Validation loss")
  plt.title('loss')
  plt.legend(loc="lower right")
  plt.savefig(directory + '/' + "Training-Validation Loss " + model_net + " " + job + " " + "batchsize " + str(batch_size) + " epoch " + str(epoch) + '.png')
  plt.clf(); plt.close()

def plot_accuracy(history, dataset, model_net, job, batch_size, epoch):
  directory= 'CGR-pcmer ' + dataset.split('/')[1] + ' Results/'
  if not os.path.exists(directory):
    os.makedirs(directory)

  plt.figure(figsize=(10,6))
  plt.plot(history.epoch,history.history['accuracy'], label = "Training accuracy")
  plt.plot(history.epoch,history.history['val_accuracy'], label = "Validation accuracy")
  plt.title('accuracy')
  plt.legend(loc="lower right")
  plt.savefig(directory + '/' + "Training-Validation Accuracy " + model_net + " " + job  + " " + "batchsize " + str(batch_size) + " epoch " + str(epoch) + '.png')
  plt.clf(); plt.close()
  
def metrics(X_test, y_test, model_net):
  y_predict = model_net.predict(X_test)
  y_maxPredict = np.arange(len(y_test))

  index = 0
  for a in y_predict:
    y_max = np.argmax(a)
    np.put(y_maxPredict,[index],[y_max])
    index += 1


  # check results
  return confusion_matrix(y_test, y_maxPredict), classification_report(y_test, y_maxPredict, digits=4)
  
  
def saveModel(dataset, model, net, type_encoder, job, batch_size, epoch):
  directory='CGR-pcmer ' + dataset.split('/')[1] + ' Models/'
  if not os.path.exists(directory):
    os.makedirs(directory)

  model_file = directory + type_encoder + "model " + net + " " + job + " " + "batchsize " + str(batch_size) + " epoch " + str(epoch) + ".keras"

  #with open(model_file, 'wb') as file:
    #pickle.dump(model, file) # dump and pickle for to store the object data to the file
    #file.close()
    
  model.save(model_file)
  print("Save Model!")
    
def plot_loss_accuracy(history, model, X_test, y_test, dataset, model_net, job, batch_size, epoch):
  plot_loss(history, dataset, model_net, job, batch_size, epoch)
  plot_accuracy(history, dataset, model_net, job, batch_size, epoch)

  scores = model.evaluate(X_test, y_test, verbose=2)
  acc = "%s: %.2f%%" % (model.metrics_names[1], scores[1]*100)
  return acc
  
def saveConfMatrixClassReport(net, test_acc, training_time, conf_matrix, class_report, dataset, type_encoder, job, val_acc, batch_size, epoch):
  directory= 'CGR-pcmer ' + dataset.split('/')[1] + ' Results/'
  if not os.path.exists(directory):
    os.makedirs(directory)

  results_model_file = directory + type_encoder + "results " + net + " " + job + " batchsize " + str(batch_size) + " epoch " + str(epoch) + ".txt"
  print("Save accuracy and classification report!")

  with open(results_model_file, 'w') as file:
    file.write('confusion matrix: \n' + str(conf_matrix) + '\n\n')
    file.write('classification report: \n' + str(class_report) + '\n')
    file.write(str(test_acc) + '\n')
    file.write(str(val_acc))
    file.write(str(training_time))
    file.close()
