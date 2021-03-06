import os
import torch
import numpy as np

import torch.nn as nn
import torchvision.transforms as transforms

from torch.utils.data import Dataset
from PIL import Image

MEAN = [0.485, 0.456, 0.406]
STD = [0.229, 0.224, 0.225]


class CNVData(Dataset):
    def __init__(self, args, mode='train'):

        'set up basic parameters for dataset'
        self.mode = mode

        'training locally or on bwccluster'
        if args.bwcluster:
            data_dir = os.path.join(args.work_dir_bwcluster, args.ws_dir)
        else:
            data_dir = args.data_dir

        'set data according to mode'
        if mode == 'train' or mode == 'inception':
            sample_id = args.sample_id
            chrom_id = args.chrom_id
        if mode == 'val' or mode == 'inception_val':
            sample_id = args.sample_id_val
            chrom_id = args.chrom_id_val
        if mode == 'test':
            sample_id = args.sample_id_test
            chrom_id = args.chrom_id_test

        'set up image path and label'
        data_dir_list = []
        helper_dir = ''

        for id in sample_id:
            for chr in chrom_id:

                if args.bwcluster:
                    helper_dir = os.path.join(data_dir, id, "chr_" + str(chr), "img_vis_50")
                else:
                    helper_dir = os.path.join(data_dir, id, "chr_" + str(chr), "img_train")

                for f in sorted(os.listdir(helper_dir)):
                    img_path = os.path.join(helper_dir, f)
                    img_label = str(f.split('_')[1])
                    img_label = int(img_label)

                    if args.exclude_no_cnvs:
                        if not img_label == 0:
                            img_tuple = (img_path, img_label)
                            data_dir_list.append(img_tuple)
                    else:
                        img_tuple = (img_path, img_label)
                        data_dir_list.append(img_tuple)

        self.data = data_dir_list

        'set up image trainsform'
        if self.mode == 'train':
            self.transform = transforms.Compose([
                #transforms.RandomHorizontalFlip(0.5),
                transforms.ToTensor(),  # (H,W,C)->(C,H,W), [0,255]->[0, 1.0] RGB->RGB
                transforms.Normalize(MEAN, STD)
            ])

        elif self.mode == 'val' or self.mode == 'test':
            self.transform = transforms.Compose([
                transforms.ToTensor(),  # (H,W,C)->(C,H,W), [0,255]->[0, 1.0] RGB->RGB
                transforms.Normalize(MEAN, STD)
            ])

        if self.mode == 'inception' or self.mode == 'inception_val':
            self.transform = transforms.Compose([
                transforms.Resize(299),
                transforms.CenterCrop(299),
                #transforms.RandomHorizontalFlip(0.5),
                transforms.ToTensor(),  # (H,W,C)->(C,H,W), [0,255]->[0, 1.0] RGB->RGB
                transforms.Normalize(MEAN, STD)
            ])

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):

        'get data'
        img_path, cls = self.data[idx]

        'read image'
        img = Image.open(img_path).convert('RGB')

        cls = np.array(cls)
        cls = torch.from_numpy(cls)
        cls = cls.long()

        return self.transform(img), cls
