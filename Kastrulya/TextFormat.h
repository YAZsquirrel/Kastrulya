#pragma once
#include <string>
#include <map>
namespace text
{	
	enum TextMods
	{
		NoMod = 0,
		Bold = 1, 
		Underline = 4,
	};
	std::map<std::string, std::string> textColors = {{"black","30"},
																	 {"red","31"},
																	 {"green","32"},
																	 {"yellow","33"},
																	 {"blue","34"},
																	 {"magenta","35"},
																	 {"cyan","36"},
																	 {"white","37"}};
	std::string colorize(const std::string& text, const std::string& color, TextMods mod, bool isForeground = false)
	{
		return "\033[" + mod + ';' + textColors[color] + 'm' + text + "\033[0m";
	}
}
