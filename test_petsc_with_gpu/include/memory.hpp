void MallocDeviceMemory(void **ptr, size_t size);

void FreeDeviceMemory(void *ptr);

void CopyToDevice(void *dst, const void *src, size_t size);