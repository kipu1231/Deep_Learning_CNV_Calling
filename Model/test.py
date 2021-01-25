import os
import torch
import parser
import models
import data
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import accuracy_score, confusion_matrix
import csv


def evaluate(model, data_loader):
    'set model to evaluate mode'
    model.eval()
    preds = []
    gts = []
    with torch.no_grad():  # do not need to calculate information for gradient during eval
        for idx, (imgs, gt) in enumerate(data_loader):

            if torch.cuda.is_available():
                imgs = imgs.cuda()
            pred = model(imgs)
            print(pred.shape)

            _, pred = torch.max(pred, dim=1)
            print(pred.shape)

            pred = pred.cpu().numpy().squeeze()
            gt = gt.numpy().squeeze()

            preds.append(pred)
            gts.append(gt)

    gts = np.concatenate(gts)
    preds = np.concatenate(preds)

    print(confusion_matrix(gts, preds))

    return accuracy_score(gts, preds)


def save_preds(model, data_loader):
    'set model to evaluate mode'
    model.eval()
    preds = []
    gts = []
    with torch.no_grad():  # do not need to calculate information for gradient during eval

        txt_file = os.path.join(args.save_dir, 'precictions.txt')
        with open(txt_file, "w+", newline="") as txtfile:
            writer = csv.writer(txtfile)

            for idx, (imgs, gt) in enumerate(data_loader):
                print(idx)
                imgs = imgs.cuda()
                pred = model(imgs)

                _, pred = torch.max(pred, dim=1)

                pred = pred.cpu().numpy().squeeze()
                gt = gt.numpy().squeeze()

                preds.append(pred)
                gts.append(gt)

                #for j in range(preds.size(0)):
                #    writer.writerow([preds[j].item()])

    gts = np.concatenate(gts)
    preds = np.concatenate(preds)

    print(gts)
    print(preds)

    return accuracy_score(gts, preds)


if __name__ == '__main__':
    args = parser.arg_parse()

    'setup GPU'
    if torch.cuda.is_available():
        torch.cuda.set_device(args.gpu)

    'prepare data_loader'
    print('===> prepare data loader ...')
    test_loader = torch.utils.data.DataLoader(Model.data.CNVData(args, mode='test'),
                                              batch_size=args.test_batch,
                                              num_workers=args.workers,
                                              shuffle=False)

    # no_cnv = 0
    # del_cnv = 0
    # dub_cnv = 0
    #
    # for idx, (imgs, gt) in enumerate(test_loader):
    #     y = gt.detach().numpy()
    #     num_zero = (y == 0).sum()
    #     num_ones = (y == 1).sum()
    #     num_twos = (y == 2).sum()
    #
    #     no_cnv = no_cnv + num_zero
    #     del_cnv = del_cnv + num_ones
    #     dub_cnv = dub_cnv + num_twos
    #
    # total = no_cnv + del_cnv + dub_cnv
    # print(total)
    # print(no_cnv/total)


    'prepare model'
    if args.model == "Net":
        model = Model.models.Net(args)
    elif args.model == "CNN_Net":
        model = Model.models.CNN_Net(args)

    model = nn.DataParallel(model)
    if torch.cuda.is_available():
        model.cuda()

    'resume save model'
    if torch.cuda.is_available():
        checkpoint = torch.load(args.resume)
    else:
        checkpoint = torch.load(args.resume, map_location=torch.device('cpu'))
    model.load_state_dict(checkpoint)

    #TODO: Write two different evaluation functions --> (1) by image (evaluate) / (2) by CNV

    acc = evaluate(model, test_loader)
    print('Testing Accuracy: {}'.format(acc))

