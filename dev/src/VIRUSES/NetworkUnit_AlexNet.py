from CCGRlib.module import *  # module.py import
from CCGRlib.functions_Net import (
    create_AlexNet, preprocessing, metrics, saveModel, 
    plot_loss_accuracy, saveConfMatrixClassReport
)

def main(args):
    # check compatibility: pcCCGR request RGB
    if args.type_encodingColour == "pcCCGR" and args.type_encoder == "Grayscale":
        raise ValueError(
            "Error: CCGR images generated with the 'pcCCGR' technique can only be trained using an 'RGB' encoder.\n"
            "-> Please correct your arguments: use --type_encoder RGB"
        )

    job = f"CCGR (k={args.kmer} T={args.threshold} {args.type_encodingColour})"

    X_data, y_data, nb_classes = preprocessing('AlexNet', args.type_encoder, args.dataset, args.kmer, job)

    batch_size = args.batch_size
    epoch = args.epochs

    X_data = np.array(X_data).reshape((-1, 227, 227, 3)).astype('float32')
    y_data = np.array(y_data)

    print('Data shape:', X_data.shape)
    print('Labels shape:', y_data.shape)
    print('Number of classes:', nb_classes)

    X_trainval, X_test, y_trainval, y_test = train_test_split(
        X_data, y_data, test_size=0.10, random_state=42, stratify=y_data
    )

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=20)
    tmp = 1
    start = time.time()
    model_AlexNet = Sequential()

    for train_index, val_index in skf.split(X_trainval, y_trainval):
        model_AlexNet = Sequential()
        X_train, X_val = X_trainval[train_index], X_trainval[val_index]
        y_train, y_val = y_trainval[train_index], y_trainval[val_index]
        print(f'Fold {tmp}:')

        with ProcessPoolExecutor(args.n_task) as e:
            e.map(create_AlexNet, range(args.n_task))
            model_AlexNet = create_AlexNet(model_AlexNet, (227, 227, 3), nb_classes)
            history = model_AlexNet.fit(
                X_train, y_train,
                batch_size=batch_size,
                epochs=epoch,
                validation_data=(X_val, y_val)
            )
            print(f'Fold {tmp} finished')
        tmp += 1

    end = time.time()
    val_acc = f"Validation accuracy {round(history.history['val_accuracy'][-1] * 100, 2)}%"
    training_time = f"Model training time of AlexNet Model with {args.type_encoder} encoder unit: {end - start:.2f} s"
    print(training_time)

    saveModel(args.dataset, model_AlexNet, 'AlexNet', args.type_encoder, job, batch_size, epoch)

    acc_test = plot_loss_accuracy(
        history, model_AlexNet, X_test, y_test, args.dataset, 'AlexNet', job, batch_size, epoch
    )

    conf_matrix, class_report = metrics(X_test, y_test, model_AlexNet)
    print('\n', conf_matrix, '\n', class_report)

    saveConfMatrixClassReport(
        'AlexNet', training_time, acc_test, conf_matrix, class_report,
        args.dataset, args.type_encoder, job, val_acc, batch_size, epoch
    )


if __name__ == "__main__":
    # --- Argparse CLI interface ---
    parser = argparse.ArgumentParser(
        description="Train AlexNet model with CCGR encoding"
    )
    parser.add_argument('--dataset', type=str, required=True, 
                        help="Dataset name (i.e., Coronaviruses, HIV1, HIV2, Dengue, HepatitisC, HepatitisB1, HepatitisB2, InfluenzaA)")
    parser.add_argument('--type_encoder', type=str, default="RGB", choices=["RGB", "Grayscale"], 
                        help="Type of encoder (default: RGB)")
    parser.add_argument('--kmer', type=int, default=4, 
                        help="Value of k for k-mer (default: 4)")
    parser.add_argument('--threshold', type=float, default=0.0, 
                        help="Threshold value (default: 0.0)")
    parser.add_argument('--type_encodingColour', type=str, default="pcCCGR", choices=["kCCGR", "pcCCGR"], 
                        help="Type of encoding colour (default: pcCCGR)")
    parser.add_argument('--batch_size', type=int, default=30, 
                        help="Batch size (default: 30)")
    parser.add_argument('--epochs', type=int, default=30, 
                        help="Number of training epochs (default: 30), to InfluenzaA use 60")
    parser.add_argument('--n_task', type=int, default=2, 
                        help="Number of tasks for parallel execution (default: 2)")
    
    args = parser.parse_args()
    main(args)