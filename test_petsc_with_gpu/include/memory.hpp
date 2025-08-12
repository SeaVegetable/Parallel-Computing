#pragma once
template <class T>
void MallocDeviceMemory(T **ptr, size_t size);

template <class T>
void FreeDeviceMemory(T *ptr);

template <class T>
void CopyToDevice(T *dst, const T *src, size_t size);