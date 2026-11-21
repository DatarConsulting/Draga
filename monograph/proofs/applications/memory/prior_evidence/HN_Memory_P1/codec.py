"""Fixed-width state codec. Only payload uint8 arrays are persistent memory."""
import itertools
import numpy as np


def candidates(budget):
    return np.array([b for b in itertools.product((2,4,6,8,10,12),repeat=4)
                     if sum(b)==budget],dtype=int)


def pack(codes,bits,budget):
    word=np.zeros(codes.shape[:2],dtype=np.uint64)
    offset=np.zeros(len(bits),dtype=np.uint64)
    for i in range(4):
        word |= codes[:,:,i].astype(np.uint64) << offset[:,None]
        offset += bits[:,i].astype(np.uint64)
    assert np.all(offset==budget)
    return np.stack([(word>>(8*j))&255 for j in range(budget//8)],axis=-1).astype(np.uint8)


def unpack(payload,bits):
    word=np.zeros(payload.shape[:2],dtype=np.uint64)
    for j in range(payload.shape[-1]):
        word |= payload[:,:,j].astype(np.uint64)<<(8*j)
    codes=np.empty((*word.shape,4),dtype=np.uint64)
    offset=np.zeros(len(bits),dtype=np.uint64)
    for i in range(4):
        mask=(np.uint64(1)<<bits[:,i].astype(np.uint64))-np.uint64(1)
        codes[:,:,i]=(word>>offset[:,None])&mask[:,None]
        offset += bits[:,i].astype(np.uint64)
    return codes
