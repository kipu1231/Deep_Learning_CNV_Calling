import os
import torch
import numpy as np

import torch.nn as nn
import torchvision.transforms as transforms

from torch.utils.data import Dataset
from PIL import Image

MEAN = [0.5, 0.5, 0.5]
STD = [0.5, 0.5, 0.5]


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
        if mode == 'train':
            sample_id = args.sample_id
            chrom_id = args.chrom_id
        if mode == 'val':
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
                helper_dir = os.path.join(data_dir, id, "chr_" + str(chr), "img_test3")

                for f in sorted(os.listdir(helper_dir)):
                    img_path = os.path.join(helper_dir, f)
                    img_label = str(f.split('_')[1])
                    img_tuple = (img_path, int(img_label))
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
