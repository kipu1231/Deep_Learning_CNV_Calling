import torch
import torch.nn as nn
from torchvision import models
import parser


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
            nn.LeakyReLU(True),
            nn.Dropout2d(0.3, True),
            nn.BatchNorm1d(512),

            nn.Linear(512, 128),
            nn.LeakyReLU(True),
            nn.Dropout2d(0.3, True),
            nn.BatchNorm1d(128),

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
            nn.ReLU(),
            nn.BatchNorm2d(32),
            #nn.MaxPool2d(kernel_size=2, stride=2),  # 64x64 -> 32x32

            # second block
            nn.Conv2d(32, 64, kernel_size=5, stride=1, padding=1),  # 131x131 -> 129x129
            nn.ReLU(),
            nn.BatchNorm2d(64),
            nn.MaxPool2d(kernel_size=2, stride=2), # 129x129 -> 64x64

            # third block
            nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),  # 64x64 -> 64x64
            nn.ReLU(),
            nn.BatchNorm2d(128),
            nn.MaxPool2d(kernel_size=2, stride=2),  # 64x64 -> 32x32

            # fourth block
            nn.Conv2d(128, 265, kernel_size=3, stride=1, padding=1),  # 32x32 -> 32x32
            nn.ReLU(),
            nn.BatchNorm2d(265),
            nn.MaxPool2d(kernel_size=2, stride=2),  # 32x32 -> 16x16

            #nn.AvgPool2d(kernel_size=4, stride=2),  # 16x16 -> 7x7

        )

        # classification
        self.classif = torch.nn.Sequential(
            nn.Linear(67840, 4096,bias=False),
            nn.ReLU(True),
            nn.BatchNorm1d(4096),
            #nn.BatchNorm1d(128),
            #nn.ReLU(True),
            nn.Linear(4096, 1000, bias=False),
            nn.ReLU(True),
            nn.BatchNorm1d(1000),

            nn.Linear(1000, 128, bias=False),
            nn.ReLU(True),
            nn.BatchNorm1d(128),

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


class Trans_Net(nn.Module):
    def __init__(self, args):
        super(Trans_Net, self).__init__()

        ''' load pretrained resnetmodel and freeze parameter '''
        model = models.resnet34(pretrained=True)
        num_features = model.fc.in_features
        model.fc = nn.Linear(num_features, 3)

        self.con_model = model

        # self.resmodel = torch.nn.Sequential(*(list(model.children())[:-2]))
        #
        # ''' declare layers used in this network'''
        # # first block
        # self.transconv1 = nn.ConvTranspose2d(512, 256, kernel_size=4, stride=2, padding=1, bias=False)
        # self.bn1 = nn.BatchNorm2d(256)
        # self.relu1 = nn.ReLU()  # 11x14 --> 22x28
        #
        # # second block
        # self.transconv2 = nn.ConvTranspose2d(256, 128, kernel_size=4, stride=2, padding=1, bias=False)
        # self.bn2 = nn.BatchNorm2d(128)
        # self.relu2 = nn.ReLU()  # 22x28 --> 44x56
        #
        # # third block
        # self.transconv3 = nn.ConvTranspose2d(128, 64, kernel_size=4, stride=2, padding=1, bias=False)
        # self.bn3 = nn.BatchNorm2d(64)
        # self.relu3 = nn.ReLU()  # 44x56 --> 88x112
        #
        # # fourth block
        # self.transconv4 = nn.ConvTranspose2d(64, 32, kernel_size=4, stride=2, padding=1, bias=False)
        # self.bn4 = nn.BatchNorm2d(32)
        # self.relu4 = nn.ReLU()  # 88x112 --> 176x224
        #
        # # fifth block
        # self.transconv5 = nn.ConvTranspose2d(32, 16, kernel_size=4, stride=2, padding=1, bias=False)
        # self.bn5 = nn.BatchNorm2d(16)
        # self.relu5 = nn.ReLU()  # 176x224 --> 352x448
        #
        # # sixth block
        # self.conv6 = nn.Conv2d(16, 9, kernel_size=1, stride=1, padding=0, bias=True)  # 352x448 --> 352x448
        #
        # self.drop = nn.Dropout2d()

    def forward(self, img):
        x = self.con_model(img)

        # x = self.relu1(self.bn1(self.transconv1(x)))
        #
        # x = self.relu2(self.bn2(self.transconv2(x)))
        #
        # x = self.relu3(self.bn3(self.transconv3(x)))
        #
        # x = self.relu4(self.bn4(self.transconv4(x)))
        #
        # x = self.relu5(self.bn5(self.transconv5(x)))
        #
        # x = self.conv6(x)
        #
        # x = self.drop(x)
        print(x)

        return x


class Deep_Variant_Net(nn.Module):
    def __init__(self, args):
        super(Deep_Variant_Net,self).__init__()

        ''' load pretrained resnetmodel and freeze parameter '''
        model_ft = models.inception_v3(pretrained=True, aux_logits=False)

        # Handle the auxilary net
        #num_ftrs = model_ft.AuxLogits.fc.in_features
        #model_ft.AuxLogits.fc = nn.Linear(num_ftrs, 3)
        # Handle the primary net
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, 3)
        nn.init.normal_(model_ft.fc.weight, mean=0, std=0.001)

        self.con_model = model_ft

    def forward(self, img):
        x = self.con_model(img)
        return x


if __name__ == '__main__':
    args = parser.arg_parse()
    model = Deep_Variant_Net(args)

    print(list(model.named_children()))