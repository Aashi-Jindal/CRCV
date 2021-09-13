import numpy as np
import struct

class save_idx:
    def __init__(self):
        pass

    def __get_dtype_index(self, d_type):

        if d_type==np.uint8:
            type_index = 0x00000800
        elif d_type==np.int8:
            type_index = 0x00000900
        elif d_type==np.uint16:
            type_index = 0x00000A00
        elif d_type==np.int16:
            type_index = 0x00000B00
        elif d_type==np.int32:
            type_index = 0x00000C00
        elif d_type==np.float32:
            type_index = 0x00000D00
        elif d_type==np.float64:
            type_index = 0x00000E00
        elif d_type==np.int64:
            type_index = 0x00000F00

        return type_index

    def __get_pack_ind(self, d_type):
        if d_type==np.uint8:
            type_index = 'B'
        elif d_type==np.int8:
            type_index = 'b'
        elif d_type==np.uint16:
            type_index = 'H'
        elif d_type==np.int16:
            type_index = 'h'
        elif d_type==np.int32:
            type_index = 'i'
        elif d_type==np.float32:
            type_index = 'f'
        elif d_type==np.float64:
            type_index = 'd'
        elif d_type==np.int64:
            type_index = 'q'

        return type_index

    def __get_magic_number(self, d_type, no_of_dims):
        magic_number = 0x00000000
        #print(d_type)
        magic_number |= self.__get_dtype_index(d_type)
        magic_number |= no_of_dims
        return magic_number

    def __get_header(self, d_type, d_shape):
        magic_number = self.__get_magic_number(d_type, len(d_shape))
        dims = list(d_shape)
        header = [magic_number] + dims
        #print(header)
        return header

    def save(self, data, fname):

        if not isinstance(data, np.ndarray):
            data = np.array(data)


        d_type = data.dtype

        header_bytes = self.__get_header(d_type,data.shape)
        pack_ind = self.__get_pack_ind(d_type)
        n_data = data.size

        with open(fname, 'wb') as fid:
            for val in header_bytes:
                fid.write(struct.pack('>I', val))
            data = data.flatten()
            for i in range(n_data):
                fid.write(struct.pack(pack_ind, data[i]))

class load_idx:
    def __init__(self):
        self.header_dtype = np.dtype(np.uint32).newbyteorder('>')

    def __get_magic_number(self):
        self.magic_number = np.frombuffer(self.fstream.read(4), dtype=self.header_dtype)
        return self.magic_number

    def __get_dtype_index(self, type_index):
        if type_index == int('0x08', 16):
            dt = np.dtype(np.uint8)
        elif type_index == int('0x09', 16):
            dt = np.dtype(np.int8)
        elif type_index == int('0x0A', 16):
            dt = np.dtype(np.uint16)
        elif type_index == int('0x0B', 16):
            dt = np.dtype(np.int16)
        elif type_index == int('0x0C', 16):
            dt = np.dtype(np.int32)
        elif type_index == int('0x0D', 16):
            dt = np.dtype(np.float32)
        elif type_index == int('0x0E', 16):
            dt = np.dtype(np.float64)
        elif type_index == int('0x0F', 16):
            dt = np.dtype(np.int64)

        return dt

    def __extract_header(self):
        mask_dim = int('0x000000ff',16)
        mask_datatype = int('0x0000ff00',16)
        no_of_dimensions = np.bitwise_and(self.magic_number, np.array(mask_dim, dtype=np.uint32))[0]
        datatype_index = np.right_shift(np.bitwise_and(self.magic_number, np.array(mask_datatype, dtype=np.uint32)),8)

        dt = self.__get_dtype_index(datatype_index)

        dimensions = np.empty(no_of_dimensions, dtype=np.uint32)

        for i in range(no_of_dimensions):
            read_val = np.frombuffer(self.fstream.read(4),dtype=self.header_dtype)
            dimensions[i] = read_val

        return dimensions, dt

    def load(self, fname):
        self.fstream = open(fname, 'rb')
        #magic_number = self.__get_magic_number()
        #print(magic_number)
        self.__get_magic_number()
        [dimensions, dt] = self.__extract_header()
        total_bytes_to_be_read = np.prod(dimensions, dtype=np.uint32)*dt.itemsize
        data = np.frombuffer(self.fstream.read(total_bytes_to_be_read),dtype=dt)
        data = np.reshape(data,dimensions)
        self.fstream.close()
        return data



if __name__ == '__main__':
    s = save_idx()
    data = np.random.random(size=(9,9))
    print(data.dtype)
    print(data.shape)
    fname = './test/data.idx'
    s.save(data, fname)

    l = load_idx()
    data_ = l.load(fname)

    print(data)
    print(data_)
