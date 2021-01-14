import vcf

all_svs = ["DEL", "INS", "DUP", "DUP:TANDEM", "INV", "ITX", "CTX"]
cnv_of_interest = ["DEL", "DUP"]


class CNV:
    def __init__(self, chrom=None, start=0, end=0, end_sup=0, name=None, sv_type=None, length=0, gt="./1",
                 info=None, sample=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.end_sup = end_sup
        self.name = name
        self.sv_type = sv_type
        self.length = length
        self.info = info
        self.gt = gt
        self.is_precise = False
        self.sample = sample

    def to_vcf_record(self, fasta_handle=None, sample=""):
        if self.start <= 0:
            return None
        if self.sv_type not in cnv_of_interest:
            return None

        # formulate the INFO field
        info = self.info
        sv_len = -self.length if self.sv_type == "DEL" else self.length
        info.update({"SVLEN": sv_len,
                     "SVTYPE": self.sv_type})
        if self.sv_type in ["DEL", "DUP", "ITX", "CTX"]:
            info["END"] = self.end

        if not self.is_precise:
            info.update({"IMPRECISE": True})

        #if self.cipos:
        #    info.update({"CIPOS": (",".join([str(a) for a in self.cipos]))})

        vcf_record = vcf.model._Record(self.chrom,
                                       self.start,
                                       ".",
                                       fasta_handle.fetch(self.chrom, max(0, self.start - 1),
                                                          max(1, self.start)) if fasta_handle else "N",
                                       [vcf.model._SV(self.sv_type)],
                                       ".",
                                       "PASS" if self.is_validated else "LowQual",
                                       info,
                                       "GT",
                                       [0],
                                       [vcf.model._Call(None, sample, vcf.model.make_calldata_tuple("GT")(GT="1/1"))])
        return vcf_record
