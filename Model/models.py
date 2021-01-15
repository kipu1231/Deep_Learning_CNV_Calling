import torch
import torch.nn as nn
from torchvision import models


class Net(nn.Module):

    def __init__(self, args):
        super(Net, self).__init__()  

        ''' declare layers used in this network'''
        self.conv_model = torch.nn.Sequential(
            nn.Conv2d(3, 96, kernel_size=11, stride=1, padding=0),  # 265x265 -> 255x255
            #nn.BatchNorm2d(32),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=3),  # 255x255 -> 85x85
            # nn.Dropout2d(0.25, True),

            nn.Conv2d(96, 256, kernel_size=5, stride=1, padding=0),  # 85x85 -> 81x81
            #n.BatchNorm2d(64),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=3),  # 81x81 -> 27x27
            #nn.Dropout2d(0.25, True),

            nn.Conv2d(256, 384, kernel_size=3, stride=1, padding=0),  # 27x27 -> 25x25
            # n.BatchNorm2d(64),
            nn.LeakyReLU(True),

            nn.Conv2d(384, 256, kernel_size=3, stride=1, padding=0),  # 25x25 -> 23x23
            #nn.BatchNorm2d(128),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=3, padding=1),  # 23x23 -> 8x8
            #nn.Dropout2d(0.3, True),

        )

        self.classifier = torch.nn.Sequential(
            nn.Linear(16384, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(True),
            nn.Dropout2d(0.3, True),

            nn.Linear(512, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(True),
            nn.Dropout2d(0.3, True),

            nn.Linear(128, 3)
        )

    def forward(self, img):
        x = self.conv_model(img)
        # nach dem Conv hat der Tensor die Größe (batchsize, channels, x, y)
        # um das ganze dem Classifier zu übergeben, reshaped man den Tensor und
        # macht exakt batchsize columns und entsprechend viele rows
        x = x.view(x.size(0), -1)  # flatten
        img_label = self.classifier(x)

        return img_label


class CNN_Net(nn.Module):

    def __init__(self, args):
        super(CNN_Net, self).__init__()

        ''' declare layers used in this network'''
        self.con_model = torch.nn.Sequential(
            # 3x265x265 --> 32x246x246

            # first block
            nn.Conv2d(3, 32, kernel_size=7, stride=2, padding=1),  # 265x265 -> 131x131
            nn.BatchNorm2d(32),
            nn.ReLU(),
            #nn.MaxPool2d(kernel_size=2, stride=2),  # 64x64 -> 32x32

            # second block
            nn.Conv2d(32, 64, kernel_size=5, stride=1, padding=1),  # 131x131 -> 129x129
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2), # 129x129 -> 64x64

            # third block
            nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),  # 64x64 -> 64x64
            nn.BatchNorm2d(128),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),  # 64x64 -> 32x32

            # fourth block
            nn.Conv2d(128, 265, kernel_size=3, stride=1, padding=1),  # 32x32 -> 32x32
            nn.BatchNorm2d(265),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),  # 32x32 -> 16x16

            #nn.AvgPool2d(kernel_size=4, stride=2),  # 16x16 -> 7x7

        )

        # classification
        self.classif = torch.nn.Sequential(
            nn.Linear(67840, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(True),

            nn.Linear(128, 3)
        )

    def forward(self, img):
        x = self.con_model(img)
        # nach dem Conv hat der Tensor die Größe (batchsize, channels, x, y)
        # um das ganze dem Classifier zu übergeben, reshaped man den Tensor und
        # macht exakt batchsize columns und entsprechend viele rows
        x = x.view(x.size(0), -1)
        img_label = self.classif(x)

        return img_label