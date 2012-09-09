#include "tcpsocket.h"
#include <algorithm>

Socket::Socket() {
}

Socket::Socket(const char* host, int port) {
	IPaddress ip;
	if(SDLNet_ResolveHost(&ip, host, port) == -1) {
		throw NetworkException();
	}
	
	if ((socket = SDLNet_TCP_Open(&ip)) == 0) {
		throw NetworkException();
	}

	set.addSocket(*this);
}

Socket::Socket(const TCPsocket& socket)
	: socket(socket) {
	set.addSocket(*this);
}

bool Socket::read(void* buffer, int length, Uint32 timeout, bool fineTimeoutControl) {
	Uint32 start = SDL_GetTicks();
	Uint32 current;

	if (timeout == 0 || !fineTimeoutControl) {
		int totalReadedBytes = 0;

		while (totalReadedBytes < length) {
			current = SDL_GetTicks();

			if (timeout != 0 && ((current - start) >= timeout || !checkForActivity(timeout - (current - start)))) {
				return false;
			}

			int readedCount;
			if ((readedCount = SDLNet_TCP_Recv(socket,reinterpret_cast<char*>(buffer) + totalReadedBytes, length - totalReadedBytes)) <= 0) {
				throw NetworkException();
			}
			totalReadedBytes += readedCount;
		}
		return true;
	} else {
		for (int i = 0; i < length; i++) {
			current = SDL_GetTicks();

			if ((current - start) >= timeout || !checkForActivity(timeout - (current - start))) {
				return false;
			}

			if (SDLNet_TCP_Recv(socket, reinterpret_cast<char*>(buffer) + i, 1) <= 0) {
				throw NetworkException();
			}
		}

		return true;
	}
}

void Socket::write(const void* buffer, int length) {
	if (SDLNet_TCP_Send(socket, buffer, length) < length) {
		throw NetworkException();
	}
}

bool Socket::readLine(char* buffer, int maxSize, Uint32 timeout) {
	if ((int) readLineBuffer.size() >= maxSize) {
		return false;
	}

	Uint32 start = SDL_GetTicks();
	Uint32 current;

	while ((int) readLineBuffer.size() < maxSize) {
		current = SDL_GetTicks();

		if (timeout != 0 && ((current - start) >= timeout || !checkForActivity(timeout - (current - start)))) {
			return false;
		}

		char ch;
		if (SDLNet_TCP_Recv(socket, &ch, 1) <= 0) {
			throw NetworkException();
		}
		if (ch == '\n') {
			strcpy(buffer, readLineBuffer.c_str());
			readLineBuffer.clear();
			break;
		} else {
			readLineBuffer += ch;
		}
	}

	return true;
}

void Socket::writeLine(const char* buffer) {
	int length = strlen(buffer);

	if (SDLNet_TCP_Send(socket, buffer, length) < length) {
		throw NetworkException();
	}
	if (SDLNet_TCP_Send(socket, "\n", 1) < 1) {
		throw NetworkException();
	}
}

void Socket::close() {
	SDLNet_TCP_Close(socket);
}

bool Socket::checkForActivity(Uint32 timeout) {
	return set.hasActiveSockets(timeout);
}

TCPsocket Socket::getSDLSocket() const {
	return socket;
}

std::string Socket::getSocketAddressAsText() {
	IPaddress* address = SDLNet_TCP_GetPeerAddress(socket);

	Uint8* octets = reinterpret_cast<Uint8*>(&address->host);
	
	char resultString[24];
	sprintf(resultString, "%d.%d.%d.%d:%d", octets[0], octets[1], octets[2], octets[3], address->port);

	return resultString;
}

bool Socket::operator==(const Socket& s) {
	return socket == s.socket;
}

ServerSocket::ServerSocket(int port) {
	IPaddress ip;
	if(SDLNet_ResolveHost(&ip, 0, port) == -1) {
		throw NetworkException();
	}
	
	if ((socket = SDLNet_TCP_Open(&ip)) == 0) {
		throw NetworkException();
	}
}

bool ServerSocket::accept(Socket& client) {
	TCPsocket clientSocket = SDLNet_TCP_Accept(socket);
	client = Socket(clientSocket);

	return clientSocket != 0;
}

void ServerSocket::close() {
	SDLNet_TCP_Close(socket);
}

SocketSet::SocketSet()
	: set(0) {
}

SocketSet::~SocketSet() {
	if (set != 0) {
		SDLNet_FreeSocketSet(set);
	}
}

void SocketSet::addSocket(const Socket& socket) {
	sockets.push_back(socket);
}

void SocketSet::removeSocket(const Socket& socket) {
	std::vector<Socket>::iterator it = std::find(sockets.begin(), sockets.end(), socket);
	if (it != sockets.end()) {
		sockets.erase(it);
	}
}

void SocketSet::prepareSet() {
	if (set != 0) {
		SDLNet_FreeSocketSet(set);
	}

	set = SDLNet_AllocSocketSet(sockets.size());
	if (!set) {
		throw NetworkException();
	}

	for (size_t i = 0; i < sockets.size(); i++) {
		if (SDLNet_TCP_AddSocket(set, sockets[i].getSDLSocket()) == -1) {
			throw NetworkException();
		}
	}
}

bool SocketSet::hasActiveSockets(Uint32 timeout) {
	prepareSet();

	int numberOfReadySockets = SDLNet_CheckSockets(set, timeout);
	if (numberOfReadySockets == -1) {
		throw NetworkException();
	}

	return numberOfReadySockets > 0;
}

std::vector<Socket> SocketSet::getActiveSockets(Uint32 timeout) {
	prepareSet();

	int numberOfReadySockets = SDLNet_CheckSockets(set, timeout);
	if (numberOfReadySockets == -1) {
		throw NetworkException();
	}

	std::vector<Socket> result;

	for (size_t i = 0; i < sockets.size() && numberOfReadySockets > 0; i++) {
		if (SDLNet_SocketReady(sockets[i].getSDLSocket())) {
			result.push_back(sockets[i]);
			numberOfReadySockets--;
		}
	}

	return result;
}
