#ifndef __TCPSOCKET_H__
#define __TCPSOCKET_H__

#include "SDL/SDL_net.h"
#include <exception>
#include <vector>
#include <string>

class NetworkException : public std::exception {
};

class Socket;

class SocketSet {
	SDLNet_SocketSet set;
	std::vector<Socket> sockets;

	void prepareSet();
public:
	SocketSet();
	~SocketSet();
	void addSocket(const Socket& socket);
	void removeSocket(const Socket& socket);
	bool hasActiveSockets(Uint32 timeout);
	std::vector<Socket> getActiveSockets(Uint32 timeout);
};

class Socket {
	TCPsocket socket;
	SocketSet set;
	std::string readLineBuffer;
public:
	Socket();
	Socket(const char* host, int port);
	Socket(const TCPsocket& socket);
	bool read(void* buffer, int length, Uint32 timeout = 0, bool fineTimeoutControl = false);
	void write(const void* buffer, int length);
	bool Socket::readLine(char* buffer, int maxSize, Uint32 timeout = 0);
	void writeLine(const char* buffer);
	void close();
	bool checkForActivity(Uint32 timeout);
	bool operator==(const Socket& s);
	std::string getSocketAddressAsText();
	TCPsocket getSDLSocket() const;
};

class ServerSocket {
	TCPsocket socket;
public:
	ServerSocket(int port);
	bool accept(Socket& client);
	void close();
};

#endif
