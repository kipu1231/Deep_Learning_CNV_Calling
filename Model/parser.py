from __future__ import absolute_import
import argparse


def arg_parse():
    parser = argparse.ArgumentParser(description='Model for calling CNVs')

    # training, validation and test data settings
    parser.add_argument('--sample_id', type=list, default=["NA20502", "NA18525"],
                        help="list of samples that should be included for training")
    parser.add_argument('--chrom_id', type=list, default=["1", "2","3", "4", "5", "6",
                                                          "7", "8", "9", "10", "11", "12"],
                        help="List of chromosomes that should be included for training")
    parser.add_argument('--sample_id_val', type=list, default=["NA19648"],
                        help="list of samples that should be included for validation")
    parser.add_argument('--chrom_id_val', type=list, default=["1", "2","3", "4", "5", "6",
                                                                "7", "8", "9", "10", "11", "12"],
                        help="List of chromosomes that should be included for validation")
    parser.add_argument('--sample_id_test', type=list, default=["NA18511"],
                        help="list of samples that should be included for testing")
    parser.add_argument('--chrom_id_test', type=list, default=["2"],
                        help="List of chromosomes that should be included for testing")

    # training parameters
    parser.add_argument('--lr', default=0.001, type=float,
                        help="initial learning rate")
    parser.add_argument('--weight-decay', default=0.0005, type=float,
                        help="initial learning rate")
    parser.add_argument('--model', type=str, default="Net",
                        help="Name of the model used for training")
    parser.add_argument('--gpu', default=0, type=int,
                        help='In homework, please always set to 0')
    parser.add_argument('--epoch', default=10, type=int,
                        help="num of validation iterations")
    parser.add_argument('--val_epoch', default=1, type=int,
                        help="num of validation iterations")
    parser.add_argument('--train_batch', default=256, type=int,
                        help="train batch size")
    parser.add_argument('--test_batch', default=256, type=int,
                        help="test batch size")


    # Datasets parameters
    parser.add_argument('--workers', default=0, type=int,
                        help="number of data loading workers (default: 4)")
    parser.add_argument('--bwcluster', type=bool, default=True,
                        help="Define if script is run on bwcluster or not")
    parser.add_argument('--work_dir_bwcluster', type=str, default="/pfs/work7/workspace/scratch",
                        help="Path of work space directory on bwcluster")
    parser.add_argument('--ws_dir', type=str, default="db0052-ma_deepsv-0/data",
                        help="Path of data workspace directory")
    parser.add_argument('--ws_model_dir', type=str, default="db0052-ma_deepsv-0/Model",
                        help="Path of model workspace directory")
    parser.add_argument('--work_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Model",
                        help="Path of working directory")
    parser.add_argument('--data_dir', type=str, default="/Users/chiarapullem/Uni/Master/MA2/Deep_Learning_CNV_Calling/Preprocessing",
                        help="Path of data directory")

    
    # resume trained model
    parser.add_argument('--resume', type=str, default='log/Net_model_best.pth.tar',
                    help="path to the trained model")
    # others
    parser.add_argument('--save_dir', type=str, default='log')
    parser.add_argument('--random_seed', type=int, default=999)

    args = parser.parse_args()

    return args
