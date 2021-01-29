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
from ray import tune
from ray.tune import CLIReporter
from ray.tune.schedulers import ASHAScheduler


def save_model(model, save_path):
    torch.save(model.state_dict(), save_path)


def train_model(config, checkpoint_dir=None):
    ''' setup GPU '''
    if torch.cuda.is_available():
        torch.cuda.set_device(args.gpu)

    ''' setup random seed '''
    np.random.seed(args.random_seed)
    torch.manual_seed(args.random_seed)
    torch.cuda.manual_seed(args.random_seed)

    ''' load dataset and prepare data loader '''
    print('===> prepare dataloader ...')

    train_loader = torch.utils.data.DataLoader(data.CNVData(args, mode='train'),
                                               batch_size=int(config["batch_size"]),
                                               num_workers=args.workers,
                                               shuffle=True)
    val_loader = torch.utils.data.DataLoader(data.CNVData(args, mode='val'),
                                             batch_size=int(config["batch_size"]),
                                             num_workers=args.workers,
                                             shuffle=False)

    ''' load model '''
    print('===> prepare model ...')
    if args.model == "Net":
        model = models.Net(args)
    elif args.model == "CNN_Net":
        model = models.CNN_Net(args)
    elif args.model == "Trans_Net":
        model = models.Trans_Net(args)

    device = "cpu"
    if torch.cuda.is_available():
        device = "cuda:0"
        if torch.cuda.device_count() > 1:
            model = nn.DataParallel(model)
    model.to(device)

    ''' define loss '''
    # weights = [0.4, 3.0, 3.0]
    # if torch.cuda.is_available():
    #     class_weights = torch.FloatTensor(weights).cuda()
    # else:
    #     class_weights = torch.FloatTensor(weights)
    #criterion = nn.CrossEntropyLoss(weight=class_weights)
    criterion = nn.CrossEntropyLoss()

    ''' setup optimizer '''
    optimizer = torch.optim.AdamW(model.parameters(), lr=config["lr"])

    ''' save configurations '''
    if checkpoint_dir:
        model_state, optimizer_state = torch.load(
            os.path.join(checkpoint_dir, "checkpoint"))
        model.load_state_dict(model_state)
        optimizer.load_state_dict(optimizer_state)

    # ''' setup tensorboard '''
    # writer = SummaryWriter(os.path.join(args.save_dir, 'train_info'))

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
            #writer.add_scalar('loss', loss.data.cpu().numpy(), iters)
            train_info += ' loss: {:.4f}'.format(loss.data.cpu().numpy())

            print(train_info)

        if epoch % args.val_epoch == 0:
            # Validation loss
            val_loss = 0.0
            val_steps = 0
            total = 0
            correct = 0
            for i, (imgs, cls) in enumerate(val_loader):
                with torch.no_grad():

                    if torch.cuda.is_available():
                        imgs, cls = imgs.cuda(), cls.cuda()

                    outputs = model(imgs)
                    _, predicted = torch.max(outputs.data, 1)
                    total += cls.size(0)
                    correct += (predicted == cls).sum().item()

                    loss = criterion(outputs, cls)
                    val_loss += loss.cpu().numpy()
                    val_steps += 1

            tune.report(loss=(val_loss / val_steps), accuracy=correct / total)

            # ''' evaluate the model '''
            # acc = test.evaluate(model, val_loader)
            # #writer.add_scalar('val_acc', acc, iters)
            #
            # # fig = plt.figure()
            # # plt.imshow(conf_mat)
            # # writer.add_figure('confusion_mat', fig)
            #
            # print('Epoch: [{}] ACC:{}'.format(epoch, acc))
            # #tune.report(accuracy=acc)


        #     ''' save best model '''
        #     if acc > best_acc:
        #         save_model(model, os.path.join(args.save_dir, 'model_best.pth.tar'))
        #         best_acc = acc
        #
        # ''' save model '''
        # save_model(model, os.path.join(args.save_dir, 'model_{}.pth.tar'.format(epoch)))


def main(num_samples=10, max_num_epochs=10, gpus_per_trial=2):

    # define search space
    config = {
        "lr": tune.loguniform(1e-6, 1e-3),
        "batch_size": tune.choice([32, 64, 128])
    }

    scheduler = ASHAScheduler(
        metric="loss",
        mode="min",
        max_t=max_num_epochs,
        grace_period=1,
        reduction_factor=2)

    reporter = CLIReporter(
        # parameter_columns=["l1", "l2", "lr", "batch_size"],
        metric_columns=["loss", "accuracy", "training_iteration"])

    result = tune.run(
        train_model,
        resources_per_trial={"gpu": gpus_per_trial},
        config=config,
        num_samples=num_samples,
        scheduler=scheduler,
        progress_reporter=reporter)

    best_trial = result.get_best_trial("loss", "min", "last")
    print("Best trial config: {}".format(best_trial.config))
    print("Best trial final validation loss: {}".format(
        best_trial.last_result["loss"]))
    print("Best trial final validation accuracy: {}".format(
        best_trial.last_result["accuracy"]))

    # if args.model == "Net":
    #     best_trained_model = models.Net(args)
    # elif args.model == "CNN_Net":
    #     best_trained_model = models.CNN_Net(args)
    # elif args.model == "Trans_Net":
    #     best_trained_model = models.Trans_Net(args)
    #
    # device = "cpu"
    # if torch.cuda.is_available():
    #     device = "cuda:0"
    #     if gpus_per_trial > 1:
    #         best_trained_model = nn.DataParallel(best_trained_model)
    # best_trained_model.to(device)
    #
    # best_checkpoint_dir = best_trial.checkpoint.value
    # model_state, optimizer_state = torch.load(os.path.join(
    #     best_checkpoint_dir, "checkpoint"))
    # best_trained_model.load_state_dict(model_state)
    #
    # test_acc = test.evaluate(best_trained_model, device)
    # print("Best trial test set accuracy: {}".format(test_acc))


if __name__ == '__main__':
    args = parser.arg_parse()

    main(num_samples=5, max_num_epochs=5, gpus_per_trial=1)
