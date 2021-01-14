import torch
import torch.nn as nn


class Net(nn.Module):

    def __init__(self, args):
        super(Net, self).__init__()  

        ''' declare layers used in this network'''
        self.conv_model = torch.nn.Sequential(
            # 3x256x256 --> 96x246x246
            nn.Conv2d(3, 96, kernel_size=11, stride=4, padding=1), #3x256x256 --> 96x246x246
            #nn.BatchNorm2d(32),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=2), #96x246x246 -->  96x82x82
            # nn.Dropout2d(0.25, True),

            # 96x82x82 --> 256x78x78
            nn.Conv2d(96, 256, kernel_size=5, stride=1, padding=1),
            #n.BatchNorm2d(64),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=2), #256x78x78 -->  256x26x26
            #nn.Dropout2d(0.25, True),

            #256x26x26 --> 384x24x24
            nn.Conv2d(256, 384, kernel_size=3, stride=1, padding=1),
            # n.BatchNorm2d(64),
            nn.LeakyReLU(True),

            #384x24x24 --> 256x22x22
            nn.Conv2d(384, 256, kernel_size=3, stride=1, padding=1),
            #nn.BatchNorm2d(128),
            nn.LeakyReLU(True),
            nn.MaxPool2d(kernel_size=3, stride=2), #256x22x22 -->  256x8x8
            #nn.Dropout2d(0.3, True),

        )

        self.classifier = torch.nn.Sequential(
            nn.Linear(9216, 512),
            nn.BatchNorm1d(512),
            nn.LeakyReLU(True),

            nn.Linear(512, 3)
        )

    def forward(self, img):
        x = self.conv_model(img)
        # nach dem Conv hat der Tensor die Größe (batchsize, channels, x, y)
        # um das ganze dem Classifier zu übergeben, reshaped man den Tensor und
        # macht exakt batchsize columns und entsprechend viele rows
        x = x.view(x.size(0), -1)
        img_label = self.classifier(x)

        return x
