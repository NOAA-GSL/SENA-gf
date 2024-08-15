import numpy as np

class MT19937:
    def __init__(self, seed=5489):
        self.nn = 312
        self.mm = 156
        self.seed_def = 5489
        self.matrix_a = -5403634167711393303
        self.um = -2147483648  # most significant 33 bits
        self.lm = 2147483647   # least significant 31 bits

        self.pi253_1 = 1.0 / (2.0**53 - 1.0)
        self.pi253 = 1.0 / (2.0**53)
        self.pi252 = 1.0 / (2.0**52)
        self.facJFM = 10.0

        self.mt = np.zeros(self.nn, dtype=np.int64)  # array for the state vector
        self.mti = self.nn + 1  # mti==nn+1 means mt[nn] is not initialized

        self.mask_bits1=0x7FFFFFFFFFFFFFFF
        self.mask_bits11=0x001FFFFFFFFFFFFF
        self.mask_bits12=0x000FFFFFFFFFFFFF
        self.mask_bits29=0x00000007FFFFFFFF
        self.mask_bits43=0x00000000001FFFFF
        self.mask_bits62=0x0000000000000003
        
    def init_genrand64(self, seed):
        self.mt[0] = seed
        for i in range(1, self.nn):
            self.mt[i] = 6364136223846793005 * (self.mt[i - 1] ^ ((self.mt[i - 1] >> 62) & self.mask_bits62)) + i
            #print(i,self.mt[i])
        self.mti = self.nn


    def init_by_array64(self, init_key):
        c1 = 3935559000370003845
        c2 = 2862933555777941757
        key_length = len(init_key)
        i = 0
        j = 0
        k = max(self.nn, key_length)

        self.init_genrand64(19650218)
        for kk in range(k):
            self.mt[i] = (self.mt[i] ^ (c1 * ((self.mt[i] >> 62) & self.mask_bits62) ^ self.mt[i])) + init_key[j] + j
            i += 1
            j += 1
            if i >= self.nn:
                self.mt[0] = self.mt[self.nn - 1]
                i = 0
            if j >= key_length:
                j = 0

        for kk in range(self.nn - 1):
            self.mt[i] = (self.mt[i] ^ (c2 * ((self.mt[i] >> 62) & self.mask_bits62) ^ self.mt[i])) - i
            i += 1
            if i >= self.nn:
                self.mt[0] = self.mt[self.nn - 1]
                i = 0

        self.mt[0] = (1 << 63)  # MSB is 1; assuring non-zero initial array


    def genrand64_int64(self):
        mag01 = [0, self.matrix_a]
        x = 0

        if self.mti >= self.nn:  # generate nn words at one time
            if self.mti == self.nn + 1:  # if init_genrand64() has not been called
                self.init_genrand64(self.seed_def)

            for i in range(self.nn - self.mm):
                x = (self.mt[i] & self.um) | (self.mt[i + 1] & self.lm)
                self.mt[i] = (self.mt[i + self.mm] ^ ((x >> 1) & self.mask_bits1) ^ mag01[x & 1])

            for i in range(self.nn - self.mm, self.nn - 1):
                x = (self.mt[i] & self.um) | (self.mt[i + 1] & self.lm)
                self.mt[i] = (self.mt[i + self.mm - self.nn] ^ ((x >> 1) & self.mask_bits1) ^ mag01[x & 1])

            x = (self.mt[self.nn - 1] & self.um) | (self.mt[0] & self.lm)
            self.mt[self.nn - 1] = (self.mt[self.mm - 1] ^ ((x >> 1) & self.mask_bits1) ^ mag01[x & 1])

            self.mti = 0

        self.mti += 1
        x = self.mt[self.mti - 1]

        x ^= ((x >> 29) & self.mask_bits29) & 6148914691236517205
        x ^= (x << 17) & 8202884508482404352
        x ^= (x << 37) & -2270628950310912
        x ^= ((x >> 43) & self.mask_bits43)

        return x


    def genrand64_real1(self):
        return ((self.genrand64_int64() >> 11) & self.mask_bits11) * self.pi253_1

    def genrand64_real2(self):
        return ((self.genrand64_int64() >> 11) & self.mask_bits11) * self.pi253

    def genrand64_real3(self):
        return (((self.genrand64_int64() >> 12) & self.mask_bits12) + 0.5) * self.pi252

    def mt19937_real1d(self, data):
        for i in range(len(data)):
            data[i] = self.genrand64_real3() * self.facJFM

    def mt19937_real2d(self, data):
        for j in range(data.shape[1]):
            for i in range(data.shape[0]):
                data[i, j] = self.genrand64_real3() * self.facJFM

    def mt19937_real3d(self, data):
        for k in range(data.shape[2]):
            for j in range(data.shape[1]):
                for i in range(data.shape[0]):
                    data[i, j, k] = self.genrand64_real3() * self.facJFM

    def mt19937_real4d(self, data):
        for m in range(data.shape[3]):
            for k in range(data.shape[2]):
                for j in range(data.shape[1]):
                    for i in range(data.shape[0]):
                        data[i, j, k, m] = self.genrand64_real3() * self.facJFM

