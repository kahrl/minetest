/*
Minetest
Copyright (C) 2015 est31 <MTest31@outlook.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef CHAT_INTERFACE_H
#define CHAT_INTERFACE_H

#include "util/container.h"
#include <string>
#include <queue>
#include "irrlichttypes.h"

enum ChatEventType {
	CET_CHAT,
	//CET_LOG,
	CET_NICK_ADD,
	CET_NICK_REMOVE,
	CET_TIME_INFO,
};

struct ChatEvent {
	ChatEvent(ChatEventType a_type,
		const std::string &a_nick,
		const std::wstring &an_evt_msg) :
	type(a_type),
	nick(a_nick),
	evt_msg(an_evt_msg)
	{}
	ChatEvent(ChatEventType a_type,
		u64 a_game_time,
		u32 a_time) :
	type(a_type),
	game_time(a_game_time),
	time(a_time)
	{}

	ChatEventType type;
	std::string nick; // CET_CHAT, CET_NICK_ADD, CET_NICK_REMOVE
	std::wstring evt_msg; // CET_CHAT
	u64 game_time; // CET_TIME_INFO
	u32 time; // CET_TIME_INFO
	//std::string log_msg; // CET_LOG
};

struct ChatInterface {
	MutexedQueue<ChatEvent> command_queue; // chat backend --> server
	MutexedQueue<ChatEvent> outgoing_queue; // server --> chat backend
};

#endif
