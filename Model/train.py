import os
import torch
import models
import data
import test
from test import *
import numpy as np
import torch.nn as nn
#from torchsummary import summary
from tensorboardX import SummaryWriter
import logging
import parser
from datetime import datetime
import matplotlib.pyplot as plt


def save_model(model, save_path):
    torch.save(model.state_dict(), save_path)


if __name__ == '__main__':

    args = parser.arg_parse()
    print(args)

    '''create directory to save trained model and other info'''
    timeObj = datetime.now()
    stamp = "%d%d%d_%d:%d" % (timeObj.year, timeObj.month, timeObj.day, timeObj.hour, timeObj.minute)
    print(stamp)

    if args.bwcluster:
        save_dir = args.model + "_lr:" + str(args.lr) + "_" + stamp
        args.save_dir = os.path.join(args.work_dir_bwcluster, args.ws_model_dir, "log", save_dir)

    if not os.path.exists(args.save_dir):
        os.makedirs(args.save_dir)

    log_name = os.path.join(args.save_dir, "training.log")
    logging.basicConfig(filename=log_name, level=logging.DEBUG, format='%(message)s')

    for arg, value in vars(args).items():
        logging.info(value)

    ''' setup GPU '''
    if torch.cuda.is_available():
        torch.cuda.set_device(args.gpu)

    ''' setup random seed '''
    np.random.seed(args.random_seed)
    torch.manual_seed(args.random_seed)
    torch.cuda.manual_seed(args.random_seed)

    ''' load model '''
    print('===> prepare model ...')
    mode_train = 'train'
    mode_val = 'val'

    if args.model == "Net":
        model = models.Net(args)
    elif args.model == "CNN_Net":
        model = models.CNN_Net(args)
    elif args.model == "Trans_Net":
        model = models.Trans_Net(args)
    elif args.model == "Deep_Variant_Net":
        model = models.Deep_Variant_Net(args)
        mode_train = 'inception'
        mode_val = 'inception_val'

    # if torch.cuda.device_count() > 1:
    #     model = DistributedDataParallel(model)

    if torch.cuda.is_available():
        model.cuda() #load model to gpu

    ''' load pretrained model '''
    if args.use_pretrained_model == 'True':
        pre_mod_dir = os.path.join(args.work_dir_bwcluster, args.ws_model_dir, "log", args.dir_pretrained_model)
        if torch.cuda.is_available():
            checkpoint = torch.load(pre_mod_dir)
        else:
            checkpoint = torch.load(pre_mod_dir, map_location=torch.device('cpu'))
        model.load_state_dict(checkpoint)

    ''' load dataset and prepare data loader '''
    print('===> prepare dataloader ...')
    train_loader = torch.utils.data.DataLoader(data.CNVData(args, mode=mode_train),
                                               batch_size=args.train_batch,
                                               num_workers=args.workers,
                                               shuffle=True)
    val_loader = torch.utils.data.DataLoader(data.CNVData(args, mode=mode_val),
                                             batch_size=args.train_batch,
                                             num_workers=args.workers,
                                             shuffle=False)

    ''' define loss '''
    # weights = [0.4, 3.0, 3.0]
    # if torch.cuda.is_available():
    #     class_weights = torch.FloatTensor(weights).cuda()
    # else:
    #     class_weights = torch.FloatTensor(weights)

    #criterion = nn.CrossEntropyLoss(weight=class_weights)
    criterion = nn.CrossEntropyLoss()

    ''' setup optimizer '''
    if args.model == "Deep_Variant_Net":
        optimizer = torch.optim.SGD(model.parameters(), lr=0.0015, momentum=0.8, weight_decay=args.weight_decay)
    else:
        optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)

    ''' setup tensorboard '''
    writer = SummaryWriter(os.path.join(args.save_dir, 'train_info'))

    ''' train model '''
    print('===> start training ...')
    iters = 0
    best_acc = 0
    for epoch in range(1, args.epoch + 1):

        model.train()

        for idx, (imgs, cls) in enumerate(train_loader):
            train_info = 'Epoch: [{0}][{1}/{2}]'.format(epoch, idx + 1, len(train_loader))
            iters += 1

            ''' move data to gpu '''
            if torch.cuda.is_available():
                imgs, cls = imgs.cuda(), cls.cuda()

            ''' Zero parameter gradients '''
            optimizer.zero_grad()  # set grad of all parameters to zero

            ''' forward path '''
            output = model(imgs)

            ''' compute loss, backpropagation, update parameters '''
            loss = criterion(output, cls)  # compute loss

            loss.backward()
            optimizer.step()  # update parameters

            ''' write out information to tensorboard '''
            writer.add_scalar('loss', loss.data.cpu().numpy(), iters)
            train_info += ' loss: {:.4f}'.format(loss.data.cpu().numpy())

            print(train_info)

        if epoch % args.val_epoch == 0:
            ''' evaluate the model '''
            acc = test.evaluate(model, val_loader)
            writer.add_scalar('val_acc', acc, iters)

            # fig = plt.figure()
            # plt.imshow(conf_mat)
            # writer.add_figure('confusion_mat', fig)

            print('Epoch: [{}] ACC:{}'.format(epoch, acc))

            ''' save best model '''
            if acc > best_acc:
                save_model(model, os.path.join(args.save_dir, 'model_best.pth.tar'))
                best_acc = acc

        ''' save model '''
        save_model(model, os.path.join(args.save_dir, 'model_{}.pth.tar'.format(epoch)))
