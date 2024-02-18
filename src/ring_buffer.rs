pub struct RingBuffer<T> {
    length: usize,
    buffer: Vec<T>,
    write_pointer: usize,
    read_pointer: usize,
}

impl<T: Copy + Default> RingBuffer<T> {
    pub fn new(length: usize) -> Self {
        let mut buff = RingBuffer {
            length: length,
            buffer: Vec::<T>::new(),
            write_pointer: 0,
            read_pointer: 0
        };
        for i in 0..length {
            buff.buffer.push(T::default());
        }
        return buff;
    }

    pub fn reset(&mut self) {
        self.write_pointer = 0;
        self.read_pointer = 0;
    }

    // `put` and `peek` write/read without advancing the indices.
    pub fn put(&mut self, value: T) {
        self.buffer[self.write_pointer] = value;
    }

    pub fn peek(&self) -> T {
        return self.buffer[self.read_pointer];
    }

    pub fn get(&self, offset: usize) -> T {
        return self.buffer[offset]
    }

    // `push` and `pop` write/read and advance the indices.
    pub fn push(&mut self, value: T) {
        self.buffer[self.write_pointer] = value;
        self.write_pointer = (self.write_pointer + 1) % self.length;
    }

    pub fn pop(&mut self) -> T {
        let val = self.buffer[self.read_pointer];
        self.read_pointer = (self.read_pointer + 1) % self.length;
        val
    }

    pub fn get_read_index(&self) -> usize {
        return self.read_pointer;
    }

    pub fn set_read_index(&mut self, index: usize) {
        self.read_pointer = index % self.length;
    }

    pub fn get_write_index(&self) -> usize {
        return self.write_pointer;
    }

    pub fn set_write_index(&mut self, index: usize) {
        self.write_pointer = index % self.length;
    }

    pub fn len(&self) -> usize {
        return self.write_pointer;
    }

    pub fn capacity(&self) -> usize {
        return self.length
    }
}
